.PHONY: conda pip clean

conda:
	conda env create -f environment.yaml

update:
	conda env update --file environment.yaml --prune

pip:
	pip install git+https://github.com/mikgroup/sigpy.git@main
	git clone git@github.com:mlazaric/Chebyshev.git
	pip install opencv-python
	pip install nibabel
	pip install tk

clean:
	rm -rf __pycache__
	rm -rf .ipynb_checkpoints
	conda env remove -n decay-modeled-cs-paper
