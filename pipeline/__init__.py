from .utils import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, \
                    load_panel_from_s3, s3_put, illumina_readgroup
from .cpipeline import dedup