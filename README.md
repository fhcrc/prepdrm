# `prepdrm`

Analysis pipeline for identifying low frequency drug resistance mutations (DRMs) samples from the Partners PrEP study.

Manuscript in preparation.

## Running

The `provisioning/` directory contains [Ansible](http://www.ansible.com) scripts for installing all required dependencies and data - look there for details.

Alternatively, use one of these options:

* On AWS: the [AMI](http://en.wikipedia.org/wiki/Amazon_Machine_Image) `ami-0ccca73c` in `us-west-2` has all necessary dependencies and control data pre-loaded. Login as `ubuntu`.
* Using [Vagrant](http://www.vagrantup.com): run `vagrant up`.
* Using the [Docker](http://www.docker.io) image `cmccoy/prepdrm`:
  `docker pull cmccoy/prepdrm; docker run -i -t -P cmccoy/prepdrm /bin/bash`.

## Data

Data from control samples in FASTQ format are available at https://s3-us-west-2.amazonaws.com/prepdrm454controls/prepdrm454controls.tar.

# License

GPL v3 - see `COPYING`
