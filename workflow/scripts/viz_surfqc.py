import os
import matplotlib, matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

matplotlib.use("Agg")

# NOTE: Workbench requires real paths when called from snakemake

def create_overlay_scene(out_dir, template_scene, t1, lh_pial, lh_white, view_plane, no_slices):
	"""
	Sets up workbench scenes for qc at each view. 
	Infers path for other hemisphere
	"""
	
	with open(template_scene, "r") as f:
		scene_contents = f.read()
		
	t1_path, t1_name = os.path.split(os.path.realpath(t1))
	
	# Get some info about nifti
	t1_nii = nib.load(os.path.realpath(t1))
	t1_shape = t1_nii.get_fdata().shape
	t1_dims = t1_nii.header["pixdim"][:3]
	t1_origin = (t1_shape[0] // 2, t1_shape[1] // 2, t1_shape[2] // 2)
	
	# Grab hemisphere and infer other hemisphere
	try: 
		pial_path, lh_pial_name = os.path.split(os.path.realpath(lh_pial))
		white_path, lh_white_name = os.path.split(os.path.realpath(lh_white))
		rh_pial_name = lh_pial_name.replace("lh", "rh")
		rh_white_name = lh_white_name.replace("lh", "rh")
	except:
		raise FileNotFoundError("Hemispheric surfaces not found...")

	# Update scene
	scene_contents = scene_contents.replace("T1w.nii.gz.path", os.path.join(t1_path, t1_name))
	scene_contents = scene_contents.replace("T1w.nii.gz.name", t1_name)
	
	scene_contents = scene_contents.replace("lh.pial.surf.gii.path", os.path.join(pial_path, lh_pial_name))
	scene_contents = scene_contents.replace("lh.pial.surf.gii.name", lh_pial_name)
	scene_contents = scene_contents.replace("lh.white.surf.gii.path", os.path.join(white_path, lh_white_name))
	scene_contents = scene_contents.replace("lh.white.surf.gii.name", lh_white_name)
	
	scene_contents = scene_contents.replace("rh.pial.surf.gii.path", os.path.join(pial_path, rh_pial_name))
	scene_contents = scene_contents.replace("rh.pial.surf.gii.name", rh_pial_name)
	scene_contents = scene_contents.replace("rh.white.surf.gii.path", os.path.join(white_path, rh_white_name))
	scene_contents = scene_contents.replace("rh.white.surf.gii.name", rh_white_name)
	
	scene_contents = scene_contents.replace("VIEW_PLANE", view_plane.upper())
	# Grab slice for each view
	if view_plane.lower() == "axial":
		t1_origin = t1_origin[2]
		slice_dim = t1_dims[2]
	elif view_plane.lower() == "parasagittal":
		t1_origin = t1_origin[0]
		slice_dim = t1_dims[0]
		scene_contents = scene_contents.replace("z = ", "x =")
	else: 
		t1_origin = t1_origin[1]
		slice_dim = t1_dims[1]
		scene_contents = scene_contents.replace ("z = ", "y = ")
		
	extent = t1_origin // 2
	slices = np.linspace(t1_origin - extent, t1_origin + extent, no_slices, dtype=int)
	
	img_fname_list = []
	# Write scenes to file
	for idx, slice_idx in enumerate(slices):
		slice_diff = slice_idx - t1_origin
		coords = slice_diff * slice_dim
		scene_contents = scene_contents.replace(f'm_sliceCoordinate{view_plane.capitalize()}">0', f'm_sliceCoordinate{view_plane.capitalize()}">{coords}')
		
		if slice_diff != 0:
			text_coords = -1 * coords if view_plane.lower() == "parasaggital" else coords
			scene_contents = scene_contents.replace("0 mm", f"{text_coords:.2f} mm")
		
		scene_fname = f"{out_dir}/{view_plane.lower()}_{idx}.scene"
		with open(scene_fname, "w") as fname:
			fname.write(scene_contents)
			
		# Reset scene
		scene_contents = scene_contents.replace(f'm_sliceCoordinate{view_plane.capitalize()}">{coords}', f'm_sliceCoordinate{view_plane.capitalize()}">0')
		scene_contents = scene_contents.replace(f"{text_coords:.2f} mm", "0 mm")
			
		# Save scene as an image
		img_fname = f"{os.path.splitext(scene_fname)[0]}.png"
		img_fname_list.append(img_fname)
		wb_cmd = f"wb_command -show-scene {scene_fname} 1 {img_fname} 774 585 -use-window-size"
		os.system(wb_cmd)
		
		# Clean up files
		os.remove(scene_fname)
		
	return img_fname_list

def arrange_scene(out_img, img_list):
	no_slices = len(img_list) // 3
	
	fig, ax = plt.subplots(3, no_slices, figsize=(8, 4), gridspec_kw={"wspace": 0, "hspace": 0.})

	for idx, img_fname in enumerate(img_list):
		img = plt.imread(img_fname)
		ax[idx // 5][idx % 5].imshow(img)
		ax[idx // 5][idx % 5].axis("off")
		os.remove(img_fname)
		
	fname = f"{out_img}"
	plt.savefig(fname, dpi=300, bbox_inches="tight")


# Main
out_dir = os.path.split(snakemake.output.report)[0]

## Generate images of each slice
img_list = []
for view_plane in ["axial", "coronal", "parasagittal"]:
	img_list.extend(create_overlay_scene(out_dir, snakemake.input.scene_template, 
										 snakemake.input.t1, snakemake.input.lh_pial, 
										 snakemake.input.lh_white, view_plane, 5))

## Arrange and save as a single file
arrange_scene(snakemake.output.report, img_list)