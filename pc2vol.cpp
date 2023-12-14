#include<iostream>
#include <fstream>      // std::ifstream

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <DGtal/shapes/WindingNumbersShape.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace DGtal;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;


int main(int argc, char**argv)
{
  CLI::App app{"pc2vol"};
  
  bool visu=false;
  app.add_flag("--visu", visu, "Enable polyscope");
  std::string filename;
  app.add_option("-i,--input", filename, "Input point cloud (ascii, todo: PLY).") ->required()->check(CLI::ExistingFile);;
  bool dryrun=false;
  app.add_flag("--dry-run", dryrun, "Just load the point cloud, compute the bounding box and output the vol size (no visu).");
  std::string outputVol="output.vol";
  app.add_option("-o,--output", outputVol, "Output VOL (default: output.vol.");
  double h=1.0;
  app.add_option("--gridstep", h, "Gridstep parameter h (default: 1.0).");

  CLI11_PARSE(app, argc, argv);
  
  if (visu)
    polyscope::init();
  
  
  
  std::ifstream ifs (filename, std::ifstream::in);
  double x,y,z,nx,ny,nz;
  std::vector<Z3i::RealPoint> dpoints;
  std::vector<Z3i::RealPoint> dnormals;
  
  size_t nbPts= 0;
  Z3i::RealPoint lower(100,100,100),upper(0,0,0);
  while (ifs.good()) {
    ifs >>x>>y>>z>>nx>>ny>>nz;
    Z3i::RealPoint p(x,y,z);
    dpoints.push_back(p);
    dnormals.push_back(Z3i::RealPoint(nx,ny,nz));
    lower = lower.inf(p);
    upper = upper.sup(p);
    if (ifs.good()) ++nbPts;
  }
  ifs.close();
  
  std::cout<<"Point cloud bbox: "<<lower<<" "<<upper<<std::endl;
  RegularPointEmbedder<Z3i::Space> pointEmbedder;
  pointEmbedder.init( h );
  Z3i::Point lowerPoint = pointEmbedder.floor( lower )- Z3i::Point::diagonal(1);
  Z3i::Point upperPoint = pointEmbedder.ceil( upper ) + Z3i::Point::diagonal(1);
  Z3i::Domain domain(lowerPoint,upperPoint);
  trace.info() <<"Digital domain = "<<domain.size()<<" " <<domain<<std::endl;
  if (dryrun)
    exit(0);
  
  Eigen::MatrixXd points(nbPts,3);
  Eigen::MatrixXd normals(nbPts,3);
  for(auto i=0; i < nbPts;++i)
  {
    auto p = dpoints[i];
    auto n = dnormals[i];
    points(i,0) = p(0);
    points(i,1) = p(1);
    points(i,2) = p(2);
    normals(i,0) = n(0);
    normals(i,1) = n(1);
    normals(i,2) = n(2);
  }
  auto pc= polyscope::registerPointCloud("input boundary points", points);
  pc->addVectorQuantity("normals", normals);
  
  trace.beginBlock("WindingNumber BVH");
  //Winding number shape
  WindingNumbersShape<Z3i::Space> wnshape(points,normals);
  Eigen::VectorXd areas = Eigen::VectorXd::Ones(points.rows());
  //areas = 1.0/(double)points.rows() * areas;
  //wnshape.setPointAreas(areas);
  trace.endBlock();
  
   
  //Winding (batched)
  size_t size = domain.size();
  Eigen::MatrixXd queries(size,3);
  auto cpt=0;
  for(const auto &vox: domain)
  {
    Eigen::RowVector3<double> p(vox[0],vox[1],vox[2]);
    p *= h;
    queries.row(cpt) = p;
    ++cpt;
  }
  
  trace.info()<<"Number of queries = "<<cpt<<std::endl;
  auto orientations = wnshape.orientationBatch(queries);
  
  //Binary Predicate
  Z3i::DigitalSet voxels(domain);
  cpt=0;
  for(const auto &voxel: domain)
  {
    if (orientations[cpt]==INSIDE)
      voxels.insertNew(voxel);
    ++cpt;
  }
  trace.info() <<"Number of voxels = "<<voxels.size()<<std::endl;

  trace.beginBlock("Exporting");
  ImageContainerBySTLVector<Z3i::Domain, unsigned char> image(domain);
  for (const auto &p: voxels)
    image.setValue(p, 128);
  VolWriter< ImageContainerBySTLVector<Z3i::Domain, unsigned char> >::exportVol(outputVol, image);
  trace.endBlock();
  
  if (visu)
  {
    //Digital surface
    Z3i::KSpace kspace;
    kspace.init(lowerPoint, upperPoint, true);
    typedef Z3i::KSpace::SurfelSet SurfelSet;
    typedef SetOfSurfels< Z3i::KSpace, SurfelSet > MySetOfSurfels;
    typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
    typedef SurfelAdjacency<Z3i::KSpace::dimension> MySurfelAdjacency;
    
    MySurfelAdjacency surfAdj( true ); // interior in all directions.
    MySetOfSurfels theSetOfSurfels( kspace, surfAdj );
    Surfaces<Z3i::KSpace>::sMakeBoundary(theSetOfSurfels.surfelSet(),
                                         kspace,
                                         voxels,
                                         lowerPoint,
                                         upperPoint);
    //Polyscope visualization
    auto surfPtr = CountedPtr<DigitalSurface< MySetOfSurfels >>(new MyDigitalSurface(theSetOfSurfels));
    auto primalSurfaceReco   = SH3::makePrimalSurfaceMesh(surfPtr);
    
    std::vector<Z3i::RealPoint> positionsReco = primalSurfaceReco->positions();
    //Fixing the embedding
    std::for_each(std::begin(positionsReco), std::end(positionsReco), [&](Z3i::RealPoint &p){p=p*h;});
    
    std::vector<std::vector<SH3::SurfaceMesh::Vertex>> facesReco;
    for(auto face= 0 ; face < primalSurfaceReco->nbFaces(); ++face)
      facesReco.push_back(primalSurfaceReco->incidentVertices( face ));
    auto psMesh = polyscope::registerSurfaceMesh("Reconstruction "+std::to_string(h), positionsReco, facesReco);
    
    polyscope::show();
  }
  return EXIT_SUCCESS;
}
