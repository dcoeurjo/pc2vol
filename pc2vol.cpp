#include<iostream>
#include <fstream>      // std::ifstream
#include <unordered_set>

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


// Global vars
Z3i::Point lowerPoint,upperPoint;
std::unordered_set<Z3i::Point> voxels;

Z3i::RealPoint lower,upper;
std::vector<Z3i::RealPoint> dpoints;
std::vector<Z3i::RealPoint> dnormals;
Eigen::MatrixXd points;
Eigen::MatrixXd normals;
size_t nbPts= 0;
bool dryrun=false;
float factor=1.0;
float h=1.0;
bool skipCGALAreas=true;

void voxelize()
{
  RegularPointEmbedder<Z3i::Space> pointEmbedder;
  pointEmbedder.init( h );
  
  lowerPoint = pointEmbedder.floor( lower )- Z3i::Point::diagonal(1);
  upperPoint = pointEmbedder.ceil( upper ) + Z3i::Point::diagonal(1);
  Z3i::Domain domain(lowerPoint,upperPoint);
  trace.info() <<"Digital domain = "<<domain.size()<<" " <<domain<<std::endl;
  
  auto pc_size_X = upper[0]-lower[0]+1;
  auto pc_size_Y = upper[1]-lower[1]+1;
  auto pc_size_Z = upper[2]-lower[2]+1;
  //DBG std::cout<<pc_size_X<<std::endl;
  //DBG std::cout<<pc_size_Y<<std::endl;
  //DBG std::cout<<pc_size_Z<<std::endl;
  auto vox_lower = domain.lowerBound();
  auto vox_upper = domain.upperBound();
  auto vox_size_X = vox_upper[0]-vox_lower[0]+1;
  auto vox_size_Y = vox_upper[1]-vox_lower[1]+1;
  auto vox_size_Z = vox_upper[2]-vox_lower[2]+1;
  //DBG std::cout<<vox_size_X<<std::endl;
  //DBG std::cout<<vox_size_Y<<std::endl;
  //DBG std::cout<<vox_size_Z<<std::endl;
  std::cout<<"pc/voxels X scale factor = "<<pc_size_X/vox_size_X<<std::endl;
  std::cout<<"pc/voxels Y scale factor = "<<pc_size_Y/vox_size_Y<<std::endl;
  std::cout<<"pc/voxels Z scale factor = "<<pc_size_Z/vox_size_Z<<std::endl;
  
  if (dryrun)
    exit(0);
  
  double mindist=10000000.0;
  for(auto i=0; i < dpoints.size(); ++i)
    for(auto j=i+1; j < dpoints.size(); ++j)
      mindist = std::min(mindist, (dpoints[i]-dpoints[j]).norm());
  std::cout<< "Mindist = "<<mindist<<std::endl;
     
  trace.beginBlock("WindingNumber BVH");
  //Winding number shape
  WindingNumbersShape<Z3i::Space> wnshape(points,normals,skipCGALAreas);
  if (skipCGALAreas)
  {
    trace.info()<<"Skipping CGAL Areas"<<std::endl;
    Eigen::VectorXd areas = Eigen::VectorXd::Ones(points.rows());
    areas = h * factor * areas;
    wnshape.setPointAreas(areas);
  }
  trace.endBlock();
  
  
  //Winding (batched)
  size_t size = domain.size();
  Eigen::MatrixXd queries(size,3);
  auto cpt=0;
  for(const auto &vox: (domain))
  {
    Eigen::RowVector3<double> p(vox[0],vox[1],vox[2]);
    p *= h;
    queries.row(cpt) = p;
    ++cpt;
  }
  
  trace.info()<<"Number of queries = "<<cpt<<std::endl;
  auto orientations = wnshape.orientationBatch(queries);
  
  //Binary Predicate
  voxels.clear();
  cpt=0;
  for(const auto &voxel: domain)
  {
    if (orientations[cpt]==INSIDE)
      voxels.insert(voxel);
    ++cpt;
  }
  trace.info() <<"Number of voxels = "<<voxels.size()<<std::endl;
}

void update()
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
  
  Z3i::DigitalSet tmpVoxels(Z3i::Domain(lowerPoint,upperPoint));
  for(const auto &v: voxels)
    tmpVoxels.insertNew(v);
  
  Surfaces<Z3i::KSpace>::sMakeBoundary(theSetOfSurfels.surfelSet(),
                                       kspace,
                                       tmpVoxels,
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
  auto psMesh = polyscope::registerSurfaceMesh("Reconstruction ", positionsReco, facesReco);
}


void callback()
{
  ImGui::SliderFloat("Gridstep", &h, 0.0,10.0);
  ImGui::SliderFloat("Area scale factor", &factor, 0.0,1.0);
  if (ImGui::Button("Voxelize"))
  {
    voxelize();
    update();
  }
}

int main(int argc, char**argv)
{
  CLI::App app{"pc2vol"};
  
  bool visu=false;
  app.add_flag("--visu", visu, "Enable polyscope");
  std::string filename;
  app.add_option("-i,--input", filename, "Input point cloud (ascii, todo: PLY).") ->required()->check(CLI::ExistingFile);;
  app.add_flag("--dry-run", dryrun, "Just load the point cloud, compute the bounding box and output the vol size (no visu).");
  
  bool flip=false;
  app.add_flag("--flip-normals", flip, "Flip input normal vectors");
  app.add_flag("--skipCGALAreas", skipCGALAreas, "Skip CGAL area estimation (use constant factor instead --factor, default: true)");

  std::string outputVol="output.vol";
  app.add_option("-o,--output", outputVol, "Output VOL (default: output.vol.");
  app.add_option("--gridstep", h, "Gridstep parameter h (default: 1.0).");
  app.add_option("--area-scale-factor", factor, "Area associatd to the samples (relative to the gridstep) (default: 1.0");

  CLI11_PARSE(app, argc, argv);
  
  if (visu)
  {
    polyscope::init();
    polyscope::state::userCallback = callback;
  }
  
  std::ifstream ifs (filename, std::ifstream::in);
  double x,y,z,nx,ny,nz;
  
  lower=Z3i::RealPoint(std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
  upper=Z3i::RealPoint (std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min());
  while (ifs.good()) {
    ifs >>x>>y>>z>>nx>>ny>>nz;
    Z3i::RealPoint p(x,y,z);
    
    //Removing duplicates
    double mindist=std::numeric_limits<double>::max();
    for(auto i=0; i < dpoints.size(); ++i)
        mindist = std::min(mindist, (dpoints[i]-p).norm());
    if (mindist < 0.00001)
      continue;
    
    dpoints.push_back(p);
    dnormals.push_back(flip? -Z3i::RealPoint(nx,ny,nz):Z3i::RealPoint(nx,ny,nz));
    lower = lower.inf(p);
    upper = upper.sup(p);
    if (ifs.good()) ++nbPts;
  }
  ifs.close();
  
  std::cout<<"Point cloud "<<dpoints.size()<<" points  bbox: "<<lower<<" "<<upper<<std::endl;
  
  points = Eigen::MatrixXd(nbPts,3);
  normals = Eigen::MatrixXd (nbPts,3);
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
  
  if (!visu)
  {
    voxelize();
    
    trace.beginBlock("Exporting");
    ImageContainerBySTLVector<Z3i::Domain, unsigned char> image(Z3i::Domain(lowerPoint,upperPoint));
    for (const auto &p: voxels)
      image.setValue(p, 128);
    VolWriter< ImageContainerBySTLVector<Z3i::Domain, unsigned char> >::exportVol(outputVol, image);
    trace.endBlock();
  }
  else
  {
    
    auto pc= polyscope::registerPointCloud("input boundary points", points);
    pc->addVectorQuantity("normals", normals);
    
    polyscope::show();
  }
  return EXIT_SUCCESS;
}
