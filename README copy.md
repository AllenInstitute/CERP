# hubR
A tool to streamline the generate of UCSC browser hubs.

## Install
```
## Either clone this github and build the R package or install via:
install.packages("/allen/programs/celltypes/workgroups/hct/NelsonJ/Home/hubR/hubR_0.1.1.tar.gz") ## Could be out-of-date
```

## Example

The track and hub bucket names must be all lowercase and only contain "-".

### Please set your own AWS credentials in setup.aws()
```
library(hubR)

## Add your credentials for API access
setup.aws(access.key = "ACCESS_KEY", secret.key = "SECRET_KEY", region="us-west-2")

## Setup the awscli-2 command line tool, default is to use Nelson Johansen's install
setup.awscli(awscli.path = "/allen/programs/celltypes/workgroups/hct/NelsonJ/Home/v2/2.4.9/bin")

## Generate signatures and write the UCSC track hub! 
## data.dir is where the bigwig files are located
hubR(track.bucket = "nhp-bg-example-track2", hub.bucket = "nhp-bg-hub2",
        species = "NHP", region="BG", type="Multiome", cluster="Subclass", genome="rheMac10",
        data.dir="/allen/programs/celltypes/workgroups/hct/NelsonJ/hubR_example/",
        output.track.file="trackDB.txt", 
        email="nelson.johansen@alleninstitute.org")
```

## Installing awscli2 for centos7
```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
./aws/install -i . -b .

aws --version
```

