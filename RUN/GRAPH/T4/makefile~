MY_PATH=../../../DEV.15
SRC=$(MY_PATH)/ORBITER/SRC
SRC2=$(MY_PATH)/ORBITER2/SRC2


GRAPH_PATH=$(SRC)/APPS/GRAPH
TOOLS_PATH=$(SRC)/TOOLS

clean:
	- rm tournament_*
	

t:
	$(GRAPH_PATH)/graph.out -n 4 -v 2 \
		-tournament \
		-W \
		-draw_poset -t \
		-draw_graphs


# 4 tournaments

t_no_sk:
	$(GRAPH_PATH)/graph.out -n 4 -v 2 \
		-tournament \
		-no_superking \
		-W -draw_poset -t \
		-draw_graphs 

# 2 tournaments


poset:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file tournament_4_poset_lvl_6.layered_graph \
		-draw tournament_4_poset_lvl_6 \
		-rad 30000 \
		-as_tournament \
		-xin 3000000 \
		-yin 3000000 \
		-xout 1000000 \
		-yout 1000000 \
		-embedded \
		-text_underneath
	pdflatex tournament_4_poset_lvl_6.tex
	open tournament_4_poset_lvl_6.pdf

poset_no_sk:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file tournament_4_no_superking_poset_lvl_6.layered_graph \
		-draw tournament_4_no_superking_poset_lvl_6 \
		-rad 30000 \
		-as_tournament \
		-xin 3000000 \
		-yin 3000000 \
		-xout 1000000 \
		-yout 1000000 \
		-embedded
	pdflatex tournament_4_no_superking_poset_lvl_6.tex
	open tournament_4_no_superking_poset_lvl_6.pdf


