# PropaNet analysis

The code in this analysis has been largely adapted from [bhi-kimlab/PropaNet](https://github.com/bhi-kimlab/PropaNet).

## Usage

For setup instructions please follow [PropaNet dependencies](https://github.com/bhi-kimlab/PropaNet#dependency).

1. To run the PropaNet analysis and reproduce the results run the following code:

    ```bash
    bash running.sh
    ```

2. To then generate the plots as presented in the paper, you will need
    Python>3.6. To install the required packages: 

    ```bash
    python3 -m pip pandas numpy networkx plotly matplotlib 
    ```

    To then generate the plots, run: 

    ```bash
    cd PropaNet_geneSets
    mkdir nets_reproduced
    cp -r ../results_reproduced/TG/* nets_reproduced/
    python3 plotNetwork.py
    ```

    The plots will now be in the folder "glists_noTFShapes".