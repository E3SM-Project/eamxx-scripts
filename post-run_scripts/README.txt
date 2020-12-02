# --------------------------
Description of post-processing scripts

Regridding:
regrid_ne1024_largememory_to_256x512.sh - fastest way to regrid if you have access to large-memory nodes on cori

regrid_multi_ne1024_files_lowrezOutput.sh & regrid_ne1024_2D_to_256x512.sh - for regridding 2D files if you do not have access to the large memory nodes.


Moving files to HPSS:
moving_scream_restart_files_to_storage.sh - Moves hX files from run directory to HPSS with checksum calculations turned on, then the simulation will verify the move for you and write a file verifying the move

moving_ne1024_files_to_storage.sh - Moves restart files from run directory to HPSS with checksum calculations turned on, then the simulation will verify the move for you and write a file verifying the move


