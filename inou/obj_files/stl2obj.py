import sys
import pymeshlab as ml

# Usage: python stl2obj.py input.stl output.obj

ms = ml.MeshSet()
ms.load_new_mesh(sys.argv[1])                          # Load STL
ms.load_filter_script('meshlab_reset_origin.mlx')     # Optional filter
ms.apply_filter_script()
ms.save_current_mesh(sys.argv[2])                      # Save as OBJ

