# pc2vol

Eg:

` ./pc2vol -i ../../bunny.pts --gridstep 2 -o output.vol --visu`

options:
* `--dry-run`: uniquement la lecture du pointcloud et ça donne la taille du volume discret (pour régler le gridstep)
* `--visu`: active la visu polyscope
* `--gridstep`: ratio sur la taille du domaine discret par rapport à la bounding box du nuage de point. Ex: si > 1, ça va sous-échantillonner l'espace  