# environment name
ENV_NAME = splicing

# create the conda environment
.PHONY: environment
environment:
    conda env create -f splicing_quant.yaml

# clone necessary repositories
.PHONY: clone
clone:
    git clone https://github.com/mrunyan1/SpliSER.git -b parallel-combine
    git clone https://github.com/mrunyan1/leafcutter.git

# install the R environment using renv
.PHONY: install_renv
install_renv:
    @echo "Make sure you have activated the environment with: conda activate $(ENV_NAME)"
    Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv'); renv::restore()"

# run everything together
.PHONY: setup
setup: environment clone
    @echo "Setup complete. Now, activate the environment with: conda activate $(ENV_NAME) and then run 'make install_renv'"

# clean up repos if necessary
.PHONY: clean
clean:
    rm -rf SpliSER leafcutter renv.lock renv
