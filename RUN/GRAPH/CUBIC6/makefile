MY_PATH=../../../DEV.15/
SRC=$(MY_PATH)/ORBITER/SRC


GRAPH_PATH=$(SRC)/APPS/GRAPH
TOOLS_PATH=$(SRC)/APPS/TOOLS

clean:
	- rm graph_6_*

g6r3:
	$(GRAPH_PATH)/graph.out -n 6 -depth 9 -regular 3 -v 2  -draw_poset \
		-draw_graphs_at_level 9 -embedded
	pdflatex graph_6_r3_rep_9_0.tex
	open graph_6_r3_rep_9_0.pdf
	pdflatex graph_6_r3_rep_9_1.tex
	open graph_6_r3_rep_9_1.pdf

poset:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file graph_6_r3_poset_lvl_9.layered_graph \
		-draw graph_6_r3_poset_lvl_9 \
		-rad 6000 \
		-as_graph \
		-xin 1000000 \
		-yin 1000000 \
		-xout 1300000 \
		-yout 1000000 \
		-embedded \
		-y_scale 0.5 \
		-scale 0.3
	pdflatex graph_6_r3_poset_lvl_9.tex
	open graph_6_r3_poset_lvl_9.pdf

tree:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file graph_6_r3_tree_lvl_9.layered_graph \
		-draw graph_6_r3_tree_lvl_9 \
		-rad 6000 \
		-as_graph \
		-xin 1000000 \
		-yin 1000000 \
		-xout 1300000 \
		-yout 1000000 \
		-embedded \
		-y_scale 0.5 \
		-scale 0.3
	pdflatex graph_6_r3_tree_lvl_9.tex
	open graph_6_r3_tree_lvl_9.pdf

g6r2:
	$(GRAPH_PATH)/graph.out -n 6 -depth 6 -regular 2 -v 2  -draw_poset \
		-draw_graphs_at_level 6 -embedded
	pdflatex graph_6_r2_rep_6_0.tex
	open graph_6_r2_rep_6_0.pdf
	pdflatex graph_6_r2_rep_6_1.tex
	open graph_6_r2_rep_6_1.pdf

poset2:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file graph_6_r2_poset_lvl_6.layered_graph \
		-draw graph_6_r2_poset_lvl_6 \
		-rad 6000 \
		-as_graph \
		-xin 1000000 \
		-yin 1000000 \
		-xout 1300000 \
		-yout 1000000 \
		-embedded \
		-y_scale 0.5 \
		-scale 0.3
	pdflatex graph_6_r2_poset_lvl_6.tex
	open graph_6_r2_poset_lvl_6.pdf

tree2:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file graph_6_r2_tree_lvl_6.layered_graph \
		-draw graph_6_r2_tree_lvl_6 \
		-rad 6000 \
		-as_graph \
		-xin 1000000 \
		-yin 1000000 \
		-xout 1300000 \
		-yout 1000000 \
		-embedded \
		-y_scale 0.5 \
		-scale 0.3
	pdflatex graph_6_r2_tree_lvl_6.tex
	open graph_6_r2_tree_lvl_6.pdf

