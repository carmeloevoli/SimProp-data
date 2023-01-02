#!/bin/bash

URL=https://raw.githubusercontent.com/CRPropa/CRPropa3-data/master/tables/PPP

wget ${URL}/xs_neutron.txt -P ./tables/
wget ${URL}/xs_proton.txt -P ./tables/
