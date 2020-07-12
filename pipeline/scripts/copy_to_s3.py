import sys

from pipeline import s3_put


for path in sys.argv[1:]:
    print(path)
    sample = path.split("/")[4]
    s3_put("omdc-data", path, prefix=f"projects/head_and_neck/samples2/{sample}")
    
