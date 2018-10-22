if [ -d ~/dark/bostroem/research ] 
then
    jupyter nbconvert --to python --template=simplepython.tpl *.ipynb

    for i in *.py
        do 
            python $i
        done
else
    echo "Error: Must ssh to dark and start sshfs before running this script"
fi

