# environment name
ENV_NAME = splicing

# create and activate the conda environment
.PHONY: environment
environment:
	conda env create -f splicing_quant.yaml

# clone necessary repositories
.PHONY: clone
clone:
	git clone https://github.com/mrunyan1/SpliSER.git -b parallel-combine
	git clone https://github.com/mrunyan1/leafcutter.git

# activate the environment (note: environment activation cannot be automated)
.PHONY: activate
activate:
	@echo "To activate the environment, run: conda activate $(ENV_NAME)"

# install the R environment using renv
.PHONY: install_renv
install_renv:
	Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv'); renv::restore()"

# run everything together
.PHONY: setup
setup: environment clone install_renv
	@echo "Setup complete. Activate the environment with: conda activate $(ENV_NAME)"

# clean up repos if necessary
.PHONY: clean
clean:
	rm -rf SpliSER leafcutter renv.lock renv
