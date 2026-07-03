# CytOrigin

`CytOrigin` is a generative model based on VAE for single-cell transcriptomic data interpolation, where two normalized expression matrices from an earlier and a later time point were used as input and the expression of the middle time point can be simulated.

<img width="4233" height="1257" alt="资源 11" src="https://github.com/user-attachments/assets/8e8e5c1a-0c31-4325-b283-2ce849fb9d5d" />

# Install and Run CytOrigin
![anndata](https://img.shields.io/badge/anndata-0.11.4-blue) ![pandas](https://img.shields.io/badge/pandas-2.2.3-orange) ![numpy](https://img.shields.io/badge/numpy-2.2.4-green) ![scanpy](https://img.shields.io/badge/scanpy-1.11.1-purple) ![seaborn](https://img.shields.io/badge/seaborn-0.13.2-pink) ![torch](https://img.shields.io/badge/torch-2.6.0-red) ![matplotlib](https://img.shields.io/badge/matplotlib-3.10.1-gray)

We provided the code used to run CytOrigin on this website. To run CytOrigin, you should create a new conda environment and install all dependency packages listed in requirements.txt
### Create and activate conda environment with requirements installed
`CytOrigin` has only been tested with Python 3.13. To use CytOrigin, please create a new conda environment using the requirements.txt as the template.
```
conda create --name cytorigin --file requirements.txt
conda activate cytorigin
```
### Install PyTorch
We developped `CytOrigin` in a CUDA 12.4 environment. But you should select the version of PyTorch that is suitable to the CUDA version of your machine. You can find the appropriate version on the [PyTorch](https://github.com/ZJUFanLab/scNiche#:~:text=version%20on%20the-,PyTorch,-and%20DGL%20website) website.
### Run CytOrigin
Please check the step-by-step tutorial in the vignettes folder:

[Step 1 Setup environment](https://github.com/ldo2zju/Cytorigin/blob/main/vignettes/Step_1_Setup_environment.ipynb)

[Step 2 Setup model and helper functions](https://github.com/ldo2zju/Cytorigin/blob/main/vignettes/Step_2_Setup_model_and_define_functions.ipynb)

[Step 3 Load and preprocess](https://github.com/ldo2zju/Cytorigin/blob/main/vignettes/Step_3_Load_and_preprocess.ipynb)

[Step 4 Simulation](https://github.com/ldo2zju/Cytorigin/blob/main/vignettes/Step_4_Simulation.ipynb)

[Step 5 Correct the interpolated expressions](https://github.com/ldo2zju/Cytorigin/blob/main/vignettes/Step_5_Correct_the_simulated_expressions.ipynb)
## About
CytOrigin is developped by Haoxue Cao. Should you have any questions, please contact Haoxue Cao at [3200102931@zju.edu.cn](3200102931@zju.edu.cn).
## License
This project is covered under the **GPL-3.0 License**.
