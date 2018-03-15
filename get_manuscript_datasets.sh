if [ "$1" = "DatasetA" ]
then
wget https://datadryad.org/bitstream/handle/10255/dryad.173985/DatasetA.tgz
else
if [ "$1" = "DatasetC" ]
then
wget https://datadryad.org/bitstream/handle/10255/dryad.173989/DatasetC.tgz
else
if [ "$1" = "DatasetB" ]
then
wget https://datadryad.org/bitstream/handle/10255/dryad.173990/DatasetB.tgz
fi
fi
fi;
tar -xvzf "$1".tgz;
