#!/bin/bash
# CJL; (cjl2007@med.cornell.edu)
# FreeSurfer stuff adapted from Glasser et al., 2013 (NeuroImage)

# define 
# subdir
subdir=$1

# reformat subject directory path  
if [ "${subdir: -1}" = "/" ]; then
	subdir=${subdir%?};
fi

# define some variables 
RESOURCES="/path/to/resources" # resources folder is included in the github account 
FS="$RESOURCES/FS" # dir. in resources folder with FreeSurfer (FS) atlases 
FSL="$RESOURCES/FSL" # dir. in resources folder with FSL (FSL) atlases 
FS_dir="/path/to/your/freesurfer_sub_dir" # freesurfer directory 
T1="/path/to/your/image/T1.nii.gz" # path to T1 image that was submitted to FreeSurfer's recon-all 
T1_brain="/path/to/your/image/T1_brain.nii.gz" # path to T1 image that was submitted to recon-all; brain extracted 
T2="/path/to/your/image/T2.nii.gz" # path to T2 image; coregistered to T1
T1_mask="/path/to/your/image/T1_mask.nii.gz" # binary mask that excludes all non-brain voxels; 
ref_vol="/path/to/your/image/T1.nii.gz" # path to a "reference" epi_volume; for me this is a single-band reference image that I am co-registering all my me-epi volumes to when I correct for head motion prior to tedana; in non multi-band data this could be the average volume of your functional run
epi2T1_inv="/path/to/your/image/inv_warp.nii.gz" # inverse warp from t1 to acpc (a warp because it includes field map corrections; otherwise could be a rigid six DOF transformation)
tedana_dir="/path/to/your/tedana_dir"

# make fine-tuning sub. folders 
mkdir "$tedana_dir"/fine_tune/ 
mkdir "$tedana_dir"/fine_tune/imgs/ 
mkdir "$tedana_dir"/fine_tune/masks/ 
mask_dir="/$tedana_dir/masks"
mkdir "$subdir"/workspace/ # make a dir.for temp files 

######################################################################################################                                                                 

# find c_ras offset between FreeSurfer surface and volume and generate matrix to transform surfaces
MatrixX=$(mri_info "$FS_dir"/mri/brain.finalsurfs.mgz | grep "c_r" | cut -d "=" -f 5 | sed s/" "/""/g)
MatrixY=$(mri_info "$FS_dir"/mri/brain.finalsurfs.mgz | grep "c_a" | cut -d "=" -f 5 | sed s/" "/""/g)
MatrixZ=$(mri_info "$FS_dir"/mri/brain.finalsurfs.mgz | grep "c_s" | cut -d "=" -f 5 | sed s/" "/""/g)

# create directories 
for mesh in fs_lr_32k fs_lr_164k ; do
	mkdir "$FS_dir"/surf/native/ 
	mkdir "$FS_dir"/surf/native/"$mesh"/ 
done

# create transformation file 
echo "1 0 0 ""$MatrixX" > "$FS_dir"/mri/c_ras.mat
echo "0 1 0 ""$MatrixY" >> "$FS_dir"/mri/c_ras.mat
echo "0 0 1 ""$MatrixZ" >> "$FS_dir"/mri/c_ras.mat
echo "0 0 0 1" >> "$FS_dir"/mri/c_ras.mat

# define tissue types 
Types="ANATOMICAL@GRAY_WHITE ANATOMICAL@PIAL"

# sweep through hemispheres 
for hemisphere in rh lh ; do 

	# set a bunch of different 
	# ways of saying left and right
	if [ $hemisphere = "lh" ] ; then
		Hemisphere="L"
		Structure="CORTEX_LEFT"
	elif [ $hemisphere = "rh" ] ; then
		Hemisphere="R"
		Structure="CORTEX_RIGHT"
	fi

	# sweep through surface types 
	for surface in pial white ; do 

		# convert from freesurfer --> .surf.gii 
		mris_convert "$FS_dir"/surf/"$hemisphere"."$surface" \
		"$FS_dir"/surf/native/"$hemisphere"."$surface".native.surf.gii \
			
		# apply previously 
		# generated affine transformation
		wb_command -surface-apply-affine \
		"$FS_dir"/surf/native/"$hemisphere"."$surface".native.surf.gii \
		"$FS_dir"/mri/c_ras.mat "$FS_dir"/surf/native/"$hemisphere"."$surface".native.surf.gii \
		
	done


	# create midthickness surfaces
	wb_command -surface-average \
	"$FS_dir"/surf/native/"$hemisphere".midthickness.native.surf.gii \
	-surf "$FS_dir"/surf/native/"$hemisphere".pial.native.surf.gii -surf \
	"$FS_dir"/surf/native/"$hemisphere".white.native.surf.gii \
		

	# convert surface metrics
	for map in sulc@sulc@Sulc \
	thickness@thickness@Thickness \
	curv@curvature@Curvature ; do

		fsname=$(echo $map | cut -d "@" -f 1)
		wbname=$(echo $map | cut -d "@" -f 2)
		mapname=$(echo $map | cut -d "@" -f 3)
			
		mris_convert -c "$FS_dir"/surf/"$hemisphere"."$fsname" \
		"$FS_dir"/surf/"$hemisphere".white "$FS_dir"/surf/native/"$hemisphere"."$wbname".native.shape.gii 
		wb_command -set-structure "$FS_dir"/surf/native/"$hemisphere"."$wbname".native.shape.gii ${Structure} 
		wb_command -set-map-names "$FS_dir"/surf/native/"$hemisphere"."$wbname".native.shape.gii -map 1 "$hemisphere"_"$mapname"
		wb_command -metric-math "var * -1" "$FS_dir"/surf/native/"$hemisphere"."$wbname".native.shape.gii -var var \
		"$FS_dir"/surf/native/"$hemisphere"."$wbname".native.shape.gii 
		wb_command -metric-palette "$FS_dir"/surf/native/"$hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE \
		-pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true 
		
	done

	# perform various thickness specific operations;  
	wb_command -metric-math "abs(thickness)" "$FS_dir"/surf/native/"$hemisphere".thickness.native.shape.gii \
	-var thickness "$FS_dir"/surf/native/"$hemisphere".thickness.native.shape.gii 
	wb_command -metric-palette "$FS_dir"/surf/native/"$hemisphere".thickness.native.shape.gii MODE_AUTO_SCALE_PERCENTAGE \
	-pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false 
	wb_command -metric-math "thickness > 0" "$FS_dir"/surf/native/"$hemisphere".roi.native.shape.gii \
	-var thickness "$FS_dir"/surf/native/"$hemisphere".thickness.native.shape.gii 
	wb_command -metric-fill-holes "$FS_dir"/surf/native/"$hemisphere".midthickness.native.surf.gii \
	"$FS_dir"/surf/native/"$hemisphere".roi.native.shape.gii "$FS_dir"/surf/native/"$hemisphere".roi.native.shape.gii
	wb_command -metric-remove-islands "$FS_dir"/surf/native/"$hemisphere".midthickness.native.surf.gii \
	"$FS_dir"/surf/native/"$hemisphere".roi.native.shape.gii "$FS_dir"/surf/native/"$hemisphere".roi.native.shape.gii  
	wb_command -set-map-names "$FS_dir"/surf/native/"$hemisphere".roi.native.shape.gii -map 1 "$Hemisphere"_ROI
	wb_command -metric-dilate "$FS_dir"/surf/native/"$hemisphere".thickness.native.shape.gii \
	"$FS_dir"/surf/native/"$hemisphere".midthickness.native.surf.gii 10 \
	"$FS_dir"/surf/native/"$hemisphere".thickness.native.shape.gii -nearest 
	wb_command -metric-dilate "$FS_dir"/surf/native/"$hemisphere".curvature.native.shape.gii \
	"$FS_dir"/surf/native/"$hemisphere".midthickness.native.surf.gii 10 \
	"$FS_dir"/surf/native/"$hemisphere".curvature.native.shape.gii -nearest 

	# convert spherical surfaces 
	for sphere in sphere.reg sphere ; do
		mris_convert "$FS_dir"/surf/"$hemisphere"."$sphere" \
		"$FS_dir"/surf/native/"$hemisphere"."$sphere".native.surf.gii 
	done

	# convert native fs_L|R meshes to fs_LR mesh
	wb_command -surface-sphere-project-unproject \
	"$FS_dir"/surf/native/"$hemisphere".sphere.reg.native.surf.gii \
	"$FS"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere.164k_fs_"$Hemisphere".surf.gii \
	"$FS"/fs_"$Hemisphere"/fs_"$Hemisphere"-to-fs_LR_fsaverage."$Hemisphere"_LR.spherical_std.164k_fs_"$Hemisphere".surf.gii \
	"$FS_dir"/surf/native/"$hemisphere".sphere.reg.reg_LR.native.surf.gii 

	# generate inflated and very inflated surfaces 
	wb_command -surface-generate-inflated \
	"$FS_dir"/surf/native/"$hemisphere".midthickness.native.surf.gii \
	"$FS_dir"/surf/native/"$hemisphere".inflated.native.surf.gii \
	"$FS_dir"/surf/native/"$hemisphere".very_inflated.native.surf.gii \
	-iterations-scale 2.5 

	# native mesh --> 164k mesh
	for surface in white midthickness pial ; do
		wb_command -surface-resample \
		"$FS_dir"/surf/native/$hemisphere."$surface".native.surf.gii \
		"$FS_dir"/surf/native/"$hemisphere".sphere.reg.reg_LR.native.surf.gii \
		"$FS"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere.164k_fs_"$Hemisphere".surf.gii \
		BARYCENTRIC "$FS_dir"/surf/native/fs_lr_164k/"$hemisphere"."$surface".164k_fs_LR.surf.gii 
	done

	# 164k mesh --> 32k mesh
	for surface in white midthickness pial ; do
		wb_command -surface-resample \
		"$FS_dir"/surf/native/fs_lr_164k/"$hemisphere"."$surface".164k_fs_LR.surf.gii \
		"$FS"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere.164k_fs_"$Hemisphere".surf.gii \
		"$FS"/"$Hemisphere".sphere.32k_fs_LR.surf.gii BARYCENTRIC \
		"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere"."$surface".32k_fs_LR.surf.gii 
		wb_command -add-to-spec-file "$FS_dir"/surf/native/fs_lr_32k/wb.spec "$Structure" \
		"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere"."$surface".32k_fs_LR.surf.gii 
	done

	# generate inflated & very inflated surfaces 
	wb_command -surface-generate-inflated \
	"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere".midthickness.32k_fs_LR.surf.gii \
	"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere".inflated.32k_fs_LR.surf.gii \
	"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere".very_inflated.32k_fs_LR.surf.gii \
	-iterations-scale .45 

	# add inflated and very inflated surfaces to the wb.spec file 
	wb_command -add-to-spec-file "$FS_dir"/surf/native/fs_lr_32k/wb.spec "$Structure" \
	"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere".inflated.32k_fs_LR.surf.gii 
	wb_command -add-to-spec-file "$FS_dir"/surf/native/fs_lr_32k/wb.spec "$Structure" \
	"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere".very_inflated.32k_fs_LR.surf.gii 

 	# sweep through the surface metrics 
	for metric in curvature sulc thickness ; do

		# resample 
		# surface metrics 
		# from native to 164k mesh
		wb_command -metric-resample \
		"$FS_dir"/surf/native/"$hemisphere"."$metric".native.shape.gii \
		"$FS_dir"/surf/native/"$hemisphere".sphere.reg.reg_LR.native.surf.gii \
		"$FS"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere.164k_fs_"$Hemisphere".surf.gii \
		BARYCENTRIC "$FS_dir"/surf/native/fs_lr_164k/"$hemisphere"."$metric".164k_fs_LR.shape.gii 

		# resample 
		# surface metrics 
		# from 164k to 32k mesh
		wb_command -metric-resample \
		"$FS_dir"/surf/native/fs_lr_164k/"$hemisphere"."$metric".164k_fs_LR.shape.gii \
		"$FS"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere.164k_fs_"$Hemisphere".surf.gii \
		"$FS"/"$Hemisphere".sphere.32k_fs_LR.surf.gii BARYCENTRIC \
		"$FS_dir"/surf/native/fs_lr_32k/"$hemisphere"."$metric".32k_fs_LR.shape.gii 

	done

done

# sweep through surface metrics 
for metric in curvature sulc thickness ; do

	# create surface metric CIFTIs
	wb_command -cifti-create-dense-timeseries \
	"$FS_dir"/surf/native/fs_lr_32k/"$metric".32k_fs_LR.dtseries.nii \
	-left-metric "$FS_dir"/surf/native/fs_lr_32k/lh."$metric".32k_fs_LR.shape.gii \
	-roi-left "$FS"/L.atlasroi.32k_fs_LR.shape.gii \
	-right-metric "$FS_dir"/surf/native/fs_lr_32k/rh."$metric".32k_fs_LR.shape.gii \
	-roi-right "$FS"/R.atlasroi.32k_fs_LR.shape.gii 

	# add metrics to spec file
	wb_command -add-to-spec-file \
	"$FS_dir"/surf/native/fs_lr_32k/wb.spec "$Structure" \
	"$FS_dir"/surf/native/fs_lr_32k/"$metric".32k_fs_LR.dtseries.nii 

done

# remove temporary metric files 
rm "$FS_dir"/surf/native/fs_lr_32k/*.shape.gii* 

# break down the MSC network templates into separate metric files; left and right hemispheres; map to native volumme & add to spec file  
wb_command -cifti-separate "$RESOURCES"/network_templates.dscalar.nii COLUMN -metric CORTEX_LEFT "$subdir"/workspace/lh_msc_templates.func.gii -roi "$FS"/L.atlasroi.32k_fs_LR.shape.gii
wb_command -cifti-separate "$RESOURCES"/network_templates.dscalar.nii COLUMN -metric CORTEX_RIGHT "$subdir"/workspace/rh_msc_templates.func.gii -roi "$FS"/R.atlasroi.32k_fs_LR.shape.gii
wb_command -metric-to-volume-mapping "$subdir"/workspace/lh_msc_templates.func.gii "$FS_dir"/surf/native/fs_lr_32k/lh.midthickness.32k_fs_LR.surf.gii "$T1"\
"$subdir"/workspace/lh_msc_templates.nii -ribbon-constrained "$FS_dir"/surf/native/fs_lr_32k/lh.white.32k_fs_LR.surf.gii "$FS_dir"/surf/native/fs_lr_32k/lh.pial.32k_fs_LR.surf.gii
wb_command -metric-to-volume-mapping "$subdir"/workspace/rh_msc_templates.func.gii "$FS_dir"/surf/native/fs_lr_32k/rh.midthickness.32k_fs_LR.surf.gii "$T1"\
"$subdir"/workspace/rh_msc_templates.nii -ribbon-constrained "$FS_dir"/surf/native/fs_lr_32k/rh.white.32k_fs_LR.surf.gii "$FS_dir"/surf/native/fs_lr_32k/rh.pial.32k_fs_LR.surf.gii
fslmaths "$subdir"/workspace/lh_msc_templates.nii -add "$subdir"/workspace/rh_msc_templates.nii "$FS_dir"/surf/native/fs_lr_32k/msc_templates_native.nii.gz
wb_command -add-to-spec-file "$FS_dir"/surf/native/fs_lr_32k/wb.spec OTHER "$FS_dir"/surf/native/fs_lr_32k/msc_templates_native.nii.gz

# make dir. for tissue masks 
mkdir "$FS_dir"/mri/tissue_masks/ 

# create white matter, csf, brain stem, cerebellum, subcortical, cortical ribbon masks 
mri_binarize --i "$FS_dir"/mri/aparc+aseg.mgz --wm --o "$FS_dir"/mri/tissue_masks/white.mgz 
mri_convert -i "$FS_dir"/mri/tissue_masks/white.mgz -o "$FS_dir"/mri/tissue_masks/white.nii.gz --like "$T1"
mri_binarize --i "$FS_dir"/mri/aparc+aseg.mgz --ventricles --match 15 --o "$FS_dir"/mri/tissue_masks/ventricles.mgz 
mri_convert -i "$FS_dir"/mri/tissue_masks/ventricles.mgz -o "$FS_dir"/mri/tissue_masks/ventricles.nii.gz --like "$T1"
mri_binarize --i "$FS_dir"/mri/aparc+aseg.mgz --match 16 --o "$FS_dir"/mri/tissue_masks/brain_stem.mgz 
mri_convert -i "$FS_dir"/mri/tissue_masks/brain_stem.mgz -o "$FS_dir"/mri/tissue_masks/brain_stem.nii.gz --like "$T1"
mri_binarize --i "$FS_dir"/mri/aparc+aseg.mgz --match 47 8 --o "$FS_dir"/mri/tissue_masks/cerebellum.mgz 
mri_convert -i "$FS_dir"/mri/tissue_masks/cerebellum.mgz -o "$FS_dir"/mri/tissue_masks/cerebellum.nii.gz --like "$T1"
mri_binarize --i "$FS_dir"/mri/aparc+aseg.mgz --match 50 49 11 10 51 12 17 53 18 54 26 58 --o "$FS_dir"/mri/tissue_masks/subcortical.mgz 
mri_convert -i "$FS_dir"/mri/tissue_masks/subcortical.mgz -o "$FS_dir"/mri/tissue_masks/subcortical.nii.gz --like "$T1"
mri_convert -i "$FS_dir"/mri/lh.ribbon.mgz -o "$FS_dir"/mri/lh.ribbon.nii.gz --like "$T1"
mri_convert -i "$FS_dir"/mri/rh.ribbon.mgz -o "$FS_dir"/mri/rh.ribbon.nii.gz --like "$T1"
fslmaths "$FS_dir"/mri/lh.ribbon.nii.gz -add "$FS_dir"/mri/rh.ribbon.nii.gz -bin "$FS_dir"/mri/tissue_masks/grey.nii.gz 

# warp grey matter ribbon into SBRef volume space; no erosion applied
applywarp --interp=nn --in="$FS_dir"/mri/tissue_masks/grey.nii.gz --ref="$ref_vol" \
-w "$epi2T1_inv" --out="$mask_dir"/grey_ero_0.nii.gz

# combine cortical ribbon, cerebellum, & subcortical structures 
fslmaths "$FS_dir"/mri/tissue_masks/grey.nii.gz -add "$FS_dir"/mri/tissue_masks/subcortical.nii.gz "$mask_dir"/grey_subcort+cerebellum_ero_0.nii.gz
fslmaths "$mask_dir"/grey_subcort+cerebellum_ero_0.nii.gz -add "$FS_dir"/mri/tissue_masks/cerebellum.nii.gz "$mask_dir"/grey_subcort+cerebellum_ero_0.nii.gz

# warp all grey matter mask into SBRef volume space; no erosion applied (used as grey matter template when fine tuning tedana component selections)
applywarp --interp=nn --in="$mask_dir"/grey_subcort+cerebellum_ero_0.nii.gz --ref="$ref_vol" -w "$epi2T1_inv" --out="$mask_dir"/grey_subcort+cerebellum_ero_0.nii.gz # overwrite high-res; mask now in functional volume space

# substract a (slightly; 2.4mm; could change this depending on your resolution) eroded brain mask from full mask; generate mask edge
fslmaths "$T1_mask" -kernel sphere 2.4 -ero "$mask_dir"/T1_brain_mask_ero.nii.gz
fslmaths "$T1_mask" -sub "$mask_dir"/T1_brain_mask_ero.nii.gz "$mask_dir"/mask_edge.nii.gz
rm "$mask_dir"/T1_brain_mask_ero.nii.gz # delete; not needed moving forward

# subtract the full ventricles mask from a dilated (slightly; 2.4mm) mask
fslmaths "$FS_dir"/mri/tissue_masks/ventricles.nii.gz -kernel sphere 2.4 -dilD "$mask_dir"/ventricles_dil.nii.gz
fslmaths "$mask_dir"/ventricles_dil.nii.gz -sub "$FS_dir"/mri/tissue_masks/ventricles.nii.gz "$mask_dir"/ventricles_edge.nii.gz
rm "$mask_dir"/ventricles_dil.nii.gz # delete; not needed moving forward

# combine mask and ventricle edge masks; generate brain_edge.nii.gz (for detecting motion-related ICA components)
fslmaths "$mask_dir"/mask_edge.nii.gz -add "$mask_dir"/ventricles_edge.nii.gz -bin "$mask_dir"/brain_edge.nii.gz

# remove some intermediate files 
rm "$mask_dir"/mask_edge.nii.gz 
rm "$mask_dir"/ventricles_edge.nii.gz 

# apply inverse of epi --> T1 warp to brain_edge.nii.gz; bring mask into into SBRef volume space;
applywarp --interp=nn --in="$mask_dir"/brain_edge.nii.gz --ref="$ref_vol" \
-w "$epi2T1_inv" --out="$mask_dir"/brain_edge.nii.gz # overwrite high-res mask 

# divide T1w by T2w; (for detecting veins)
fslmaths "$T1" -div "$T2" -mas "$T1_mask" "$mask_dir"/veins.nii.gz -odt float 
fslmaths "$mask_dir"/veins.nii.gz -div `fslstats "$mask_dir"/veins.nii.gz -k "$T1_mask" -P 50` \
-mul 2.18 -thr 10 -min 50 -div 50 "$mask_dir"/veins.nii.gz
		
# apply inverse of epi --> T1 warp to veins.nii.gz & warp mask into into SBRef volume space;
applywarp --interp=nn --in="$mask_dir"/veins.nii.gz --ref="$ref_vol" -w "$epi2T1_inv" \
--out="$mask_dir"/veins.nii.gz # overwrite high-res mask 

# smooth components in the volume; only used for visualization 
fslmaths "$tedana_dir"/betas_OC.nii -s 2.55 "$tedana_dir"/fine_tune/betas_OC_s2.55.nii.gz # 
	
# warp MSC templates to subject's SBRef functional volume space 
applywarp --interp=spline --in="$FS_dir"/surf/native/fs_lr_32k/msc_templates_native.nii.gz --ref="$ref_vol" \
-w "$epi2T1_inv" --out="$tedana_dir"/fine_tune/network_templates.nii.gz

# T1_acpc --> func. space; used for overlaying ica components on anatomy (same for T1_acpc_bet as well)
applywarp --interp=nn --in="$T1" --ref="$ref_vol" -w "$epi2T1_inv" --out="$tedana_dir"/fine_tune/T1_func.nii.gz
applywarp --interp=nn --in="$T1_brain" --ref="$ref_vol" -w "$epi2T1_inv" --out="$tedana_dir"/fine_tune/T1_bet_func.nii.gz


                                                                 
