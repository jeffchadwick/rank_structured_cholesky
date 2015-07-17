
PYTHON=python
export PYTHONPATH=$(HOME)/local/rsc/lib/python

UFBASE=http://www.cise.ufl.edu/research/sparse/mat
RUNNER=$(PYTHON) runner.py
CSUB=csub -f pac.txt

PHONY: download setup run run_1138_bus run_Flan 1565 \
       sub_1138_bus sub_Flan_1565

setup: 1138_bus/1138_bus.mat Flan_1565/Flan_1565.mat

1138_bus/1138_bus.mat:
	mkdir -p 1138_bus
	( cd 1138_bus ; wget $(UFBASE)/HB/1138_bus.mat )

Flan_1565/Flan_1565.mat:
	mkdir -p Flan_1565
	( cd Flan_1565 ; wget $(UFBASE)/Janna/Flan_1565.mat )

run: run_1138_bus run_Flan_1565

run_1138_bus: 1138_bus/1138_bus.mat
	( cd 1138_bus ; $(PYTHON) ../runner.py run HB/1138_bus )

run_Flan_1565: Flan_1565/Flan_1565.mat
	( cd Flan_1565 ; $(PYTHON) ../runner.py Janna/Flan_1565 )

sub_1138_bus:
	$(CSUB) make run_1138_bus

sub_Flan_1565:
	$(CSUB) make run_Flan_1565

clean:
	rm csub*

distclean: clean
	rm -rf Flan_1565*
	rm -rf 1138_bus*
