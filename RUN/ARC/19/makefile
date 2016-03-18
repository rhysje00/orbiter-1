MY_PATH=../../../..
#MY_PATH=/scratch/betten/COMPILE/2014_11_01
SRC=$(MY_PATH)/ORBITER/SRC

HYPEROVAL_PATH=$(SRC)/APPS/HYPEROVAL
TOOLS_PATH=$(SRC)/APPS/TOOLS
PROJECTIVE_SPACE_PATH=$(SRC2)/PROJECTIVE_SPACE

ORDER=19
TARGET=20


all: arc





clean:
	- rm arc_$(ORDER)* report*  
	- rm -r ISO
	- rm -r STARTERS
	- rm -r LIFT_DATA
	- rm -r SOLUTIONS
	- rm hyperovals_*
	- rm classification*
	- rm PG_*



cheat:
	$(PROJECTIVE_SPACE_PATH)/cheat_sheet_PG.out -n 2 -q $(ORDER)
	pdflatex PG_2_$(ORDER).tex
	open PG_2_$(ORDER).pdf


arc:
	$(HYPEROVAL_PATH)/hyperoval.out -v 4 -q $(ORDER) \
		-target_size $(TARGET) \
		-compute_starter  \
		-starter_size $(TARGET) \
		-prefix_starter STARTER/ \
		-w 


#		-W -T \
#		-draw_poset



draw:
	$(TOOLS_PATH)/layered_graph_main.out -v 4 \
		-file arc_$(ORDER)_poset_lvl_$(TARGET).layered_graph \
		-draw arc_$(ORDER)_poset_lvl_$(TARGET) \
		-xin 10000 -yin 10000 \
		-xout 1000000 -yout 1000000 \
		-rad $(RAD) \
		-embedded
	mpost arc_$(ORDER)_poset_lvl_$(TARGET).mp
	#pdflatex arc_$(ORDER)_poset_lvl_$(TARGET).tex
	#open arc_$(ORDER)_poset_lvl_$(TARGET).pdf



