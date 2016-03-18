MY_PATH=../../../../../DEV.15
#MY_PATH=/scratch/betten/COMPILE/2015_06_13
SRC=$(MY_PATH)/ORBITER/SRC
SRC2=$(MY_PATH)/ORBITER2/SRC2


BLT_PATH=$(SRC)/APPS/BLT
TOOLS_PATH=$(SRC)/APPS/TOOLS


Q=13
Q_PLUS_ONE=14
DEPTH=5
DEPTH_MINUS_ONE=4


all: starter lift build_db read_solutions compute_orbits isomorph report


clean:
	- rm x.log
	- rm -rf G
	- rm BLT_$(Q)*
	- rm -rf GRAPHS ISO SOLUTIONS 
	- rm report_*
	- rm graph_bit*
	- rm graph_0.magma
	- rm decomp_*
	- rm classification*
	- rm decomposition*
	- rm graph_bitmatrix*
	- rm -rf STARTER_DIR
	- rm -rf G
	- rm Down*
	- rm Nb_*
	- rm log*


deep:
	$(BLT_PATH)/blt.out -v 5 -BLT -depth $(Q_PLUS_ONE) \
		-w -q $(Q) -starter -schreier $(Q_PLUS_ONE) \
		-draw_poset



draw:
	$(TOOLS_PATH)/layered_graph_main.out -v 2 \
		-file BLT_11_poset_lvl_$(Q_PLUS_ONE).layered_graph \
		-draw BLT_$(Q) \
		-xin 1000000 -yin 1000000 \
		-xout 1000000 -yout 1000000 \
		-y_stretch 1.1 \
		-rad 2000 \
		-nodes_empty \
		-embedded
	pdflatex BLT_$(Q).tex
	open BLT_$(Q).pdf


starter:
	- mkdir STARTER_DIR
	$(BLT_PATH)/blt.out -v 4 -depth $(DEPTH) \
		-W -BLT -q $(Q) \
		-starter -starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-schreier $(DEPTH) 


conv:
	/sw/bin/convert -density 600x600 -trim view.ps BLT_11_tree.png


graphs:
	- mkdir G
	$(BLT_PATH)/blt.out -v 2 -BLT   -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-create_graphs 0 1 $(DEPTH_MINUS_ONE) \
		-lex \
		-output_prefix G/




solve:
		-mkdir SOLUTIONS
		$(TOOLS_PATH)/all_rainbow_cliques.out -v 2 \
			-list_of_cases GRAPHS/list_of_cases_BLT_$(Q)_$(DEPTH)_0_1.txt \
				GRAPHS/graph_BLT_$(Q)_$(DEPTH)_%ld.bin \
			-output_file SOLUTIONS/solutions_$(Q)_$(DEPTH)_all.txt 
			

#-draw






lift:
	- mkdir G
	$(BLT_PATH)/blt.out -v 2 -BLT   -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-create_graphs 0 1 $(DEPTH_MINUS_ONE) \
		-lex \
		-output_prefix G/
	-mkdir SOLUTIONS
	$(TOOLS_PATH)/all_rainbow_cliques.out -v 2 \
		-list_of_cases G/list_of_cases_BLT_$(Q)_$(DEPTH)_0_1.txt \
		G/graph_BLT_$(Q)_$(DEPTH)_%ld.bin \
		-output_file SOLUTIONS/solutions_$(Q)_$(DEPTH)_0_1.txt
	- rm -rf G








build_db:
	$(BLT_PATH)/blt.out -v 1 -BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-build_db 




read_solutions:
	$(BLT_PATH)/blt.out -v 2 -BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-solution_directory_name SOLUTIONS/ \
		-read_solutions

read_solutions_50:
	$(BLT_PATH)/blt.out -v 2 -BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-solution_directory_name SOLUTIONS/ \
		-read_solutions_after_split 50



compute_orbits:
	$(BLT_PATH)/blt.out -v 1 -BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-compute_orbits 


isomorph:
	$(BLT_PATH)/blt.out -v 2 -BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-isomorph_testing  \
		-print_interval 250



report:
	- $(BLT_PATH)/blt.out -v 1 -BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-report 
	pdflatex report_BLT_$(Q).tex
	pdflatex report_BLT_$(Q).tex
	pdflatex report_BLT_$(Q).tex




# do it three times in order to get table of contents right



down:
	$(BLT_PATH)/blt.out -v 2 \
		-BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-down_orbits 


sample:
	$(BLT_PATH)/sample.out -v 1 -BLT -q $(Q) -depth $(DEPTH) -number_of_trials 10 -lex -seed 1 \
		-out_prefix_graphs ./GRAPHS/ \
		-out_prefix_solutions ./SOLUTIONS/ 


classification_graph:
	$(BLT_PATH)/blt.out -v 2 \
		-BLT -q $(Q) \
		-starter_depth $(DEPTH) \
		-starter_directory_name STARTER_DIR/ -starter_prefix BLT_$(Q) \
		-classification_graph $(DEPTH) 

draw_decomp:
	cp ISO/decomposition_matrix.csv .
	$(TOOLS_PATH)/draw_matrix.out -csv decomposition_matrix.csv -fname_out decomp
	pdflatex decomp_bitmatrix.tex
	open decomp_bitmatrix.pdf

draw_classification_graph:
	cp ISO/classification_graph.layered_graph .
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file classification_graph.layered_graph \
		-draw classification_graph \
		-rad 1000 \
		-xin 1000000 \
		-yin 1000000 \
		-xout 1000000 \
		-yout 1000000 \
		-embedded \
		-nodes_empty
	#mpost classification_graph.mp
	pdflatex classification_graph.tex
	open classification_graph.pdf


draw_graph:
	$(TOOLS_PATH)/draw_colored_graph.out -v 2 -export_magma graph_0.magma -file G/graph_BLT_$(Q)_$(DEPTH)_0.bin -out_file graph -coordinates 1000000 1000000 500000 500000
	pdflatex graph_bitmatrix.tex
	open graph_bitmatrix.pdf





#####

