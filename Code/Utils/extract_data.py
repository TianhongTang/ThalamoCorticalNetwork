from pathlib import Path
import shutil

# Path
src_dir = Path("/weka/scratch/ilyamonosov/Tianhong/Data/Working/GLM")
dst_dir = Path("/weka/scratch/ilyamonosov/Tianhong/Data/Extracted")
substring = "dur15_0_DeltaPure_L2=0_2_epoch3000"

# check path existance
dst_dir.mkdir(parents=True, exist_ok=True)

# loop over files
count = 0
for file_path in src_dir.iterdir():
    if file_path.is_file() and substring in file_path.name:
        # print(f"Copying {file_path} to {dst_dir}")
        shutil.copy2(file_path, dst_dir / file_path.name)
        count += 1

print(f"Done. Copied {count} files.")