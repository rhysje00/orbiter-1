MY_PATH=../../../DEV.15
SRC=$(MY_PATH)/ORBITER/SRC
SRC2=$(MY_PATH)/ORBITER2/SRC2


GRAPH_PATH=$(SRC)/APPS/GRAPH
TOOLS_PATH=$(SRC)/APPS/TOOLS
DRAW_PATH=$(SRC2)/DRAW

clean:
	- rm graph_*


petersen:
	$(GRAPH_PATH)/graph.out -n 10 -regular 3 -depth 15 -girth 5 -v 2 -W -t -x 1000000 -y 500000


tree_draw_circle:
	$(TOOLS_PATH)/treedraw.out  -file graph_10_r3_g5_15.tree  \
		-x 2000000 -y 2000000 \
		-circle  -rad 13000 -on_circle \
		-graph 10 -no_circletext
	mpost graph_10_r3_g5_15.mp
	latex view.tex 
	dvips view.dvi -T 50cm,50cm -o
	open view.ps
	#pdflatex graph_10_r3_g5_15.tex
	#open graph_10_r3_g5_15.pdf

tree_draw:
	$(TOOLS_PATH)/treedraw.out  -file graph_10_r3_g5_15.tree  \
		-x 2000000 -y 2000000 \
		-circle  -rad 18000  \
		-graph 10 -no_circletext
	mpost graph_10_r3_g5_15.mp
	latex view.tex 
	dvips view.dvi -T 50cm,50cm -o
	open view.ps
	#pdflatex graph_10_r3_g5_15.tex
	#open graph_10_r3_g5_15.pdf


id:
	$(GRAPH_PATH)/graph.out -n 10 -regular 3 -depth 15 -girth 5 -v 2 \
		-identify 0 3 4 9 13 17 21 24 28 34 36 37 40 41 43 -1 


tree:
	$(TREEDRAW_PATH)/treedraw.out -on_circle \
		-file graph_10_r3_g5_15.tree \
		-graph 10 \
		-rad 30000 \
		-x 3000000 \
		-y 3000000 \
		-no_circletext \
		-graph_perm 10 1 0 6 2 5 4 8 9 7 3 \
		-multiple_circles 2
	mpost graph_10_r3_g5_15.mp

draw_petersen:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file aux_poset_lvl_15.layered_graph \
		-draw petersen \
		-rad 80
		mpost petersen.mp


conv:
	/sw/bin/convert -density 600x600 -trim view.ps petersen_tree.png


