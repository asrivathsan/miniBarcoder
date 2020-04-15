eval "$(conda shell.bash hook)"
conda config --add channels bioconda;
conda create -n $1 python=2.7 mafft racon graphmap=0.5.2 blast seqtk git fasta3;
conda activate $1;
conda install -c anaconda biopython;

if [ $2 = '.' ]; then
	echo "cloning to diretory miniBarcoder"
	git clone https://github.com/asrivathsan/miniBarcoder/;
    cd miniBarcoder;
    python setup.py install; 
else
	echo "cloning to diretory $2";
	git clone https://github.com/asrivathsan/miniBarcoder/ $2;
    cd $2;
    python setup.py install; 
fi;



if command -v mafft >/dev/null 2>&1 ; then
    echo "mafft is found"
else
    echo "mafft not found"
fi;

if command -v blastn >/dev/null 2>&1 ; then
    echo "blast is found"
else
    echo "blast not found"
fi;
if command -v glsearch36 >/dev/null 2>&1 ; then
    echo "glsearch36 is found"
else
    echo "glsearch36 not found"
fi;
if command -v racon >/dev/null 2>&1 ; then
    echo "racon is found"
else
    echo "racon not found"
fi;
if command -v graphmap >/dev/null 2>&1 ; then
    echo "graphmap is found"
else
    echo "graphmap not found"
fi;
if command -v seqtk >/dev/null 2>&1 ; then
    echo "seqtk is found"
else
    echo "seqtk not found"
fi;

