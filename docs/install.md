[Back to main page](../README.md)  

# Installing Mustache
First, install miniconda3. This is an environment management system that should keep everything organized.

Once installed, clone this github directory to some location where it can be stored permanently.

    git clone https://github.com/durrantmm/mustache.git
    
Then enter the newly downlaoded mustache directory

    cd mustache
    
Install the mustache conda environment and other dependencies with

    bash install.sh

Now activate the environment with
    
    source activate mustache
    
This is a step that must be repeated whenever using mustache from within this environment.

Now install mustache from with the command

    pip install .
    
Once complete, you can check to see if mustache installed properly by simply typing

    mustache
   
This can then be called from anywhere on the file system while in the `mustache` conda environment.


[NEXT: Install or update software](tutorial.md)



