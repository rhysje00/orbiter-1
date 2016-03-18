MY_PATH=../../../DEV.15/
SRC=$(MY_PATH)/ORBITER/SRC


CODES_PATH=$(SRC)/APPS/CODES
TOOLS_PATH=$(SRC)/APPS/TOOLS


clean:
	- rm codes_*
	- rm Twisted*
	- rm log

18_9_8:
	$(CODES_PATH)/codes.out -v 2 -n 18 -k 9 -q 4 -d 8 -lex -w -draw_poset


# Jan 2, 2015 new MacBook: 4 min 48 sec
# Aug 23, 2015 new MacBook: 4 min 35 sec


18_9_7:
	$(CODES_PATH)/codes.out -v 2 -n 50 -k 41 -q 4 -d 7 -w 
 
draw:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file codes_18_9_4_8_poset_lvl_18.layered_graph \
		-draw TwistedTensor \
		-xin 10000 -yin 10000 \
		-xout 2000000 -yout 2000000 \
		-y_stretch 0.6 \
		-scale .25 \
		-line_width 1.0 \
		-nodes_empty \
		-rad 25 \
		-embedded
	pdflatex TwistedTensor.tex
	open TwistedTensor.pdf

