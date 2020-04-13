from boto3 import client
import os



def s3_list_keys(bucket, prefix, extension=""):
    """ Returns a dict of all objects in bucket that have the specified prefix and extension.
    """
    if extension:
        extension = ".{}".format(extension)
    s3 = client("s3")
    response = {}
    kwargs = {}
    keys = {}
    while response.get("IsTruncated", True):
        response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix, **kwargs)
        for content in response.get("Contents", ()):
            if content["Key"].endswith(extension):
                keys[content["Key"]] = content
        kwargs = {"ContinuationToken": response.get("NextContinuationToken", None)}
    return keys



def s3_get(bucket, key, filename):
    s3 = client("s3")
    s3.download_file(bucket, key, filename)



for key in s3_list_keys("omdc-data", "projects/accept", extension="annotation.tsv"):
    fn = key.split("/")[-1]
    if "-c-" in fn:
        print(fn)
        s3_get("omdc-data", key, os.path.join("/home/ed/Data/accept/annotation", fn))


















