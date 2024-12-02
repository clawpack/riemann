#Collect Pytest tests for AMRClaw
cd ${CLAW}/amrclaw/tests
pathlist=()
for dir in $(ls -d */); do
    cd $dir
    #Check if there are any test files
    #Capture the error message and check if it is empty
    if [ -z "$(ls *test*.py 2>/dev/null )" ]; then
        cd ..
        continue
    fi
    for file in $(ls *test*.py); do
        pathlist+=($(pwd)/$file)
    done
    cd ..
done

#Run the tests
pytest ${pathlist[@]}