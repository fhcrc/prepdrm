#!/bin/sh

set -e
set -u

sudo apt-get update -qq
sudo apt-get install -y -qq python-dev python-pip python-apt python-pycurl
sudo pip install ansible==1.5.4
