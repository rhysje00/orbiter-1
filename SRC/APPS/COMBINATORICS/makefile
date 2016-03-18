SRC=../../../SRC
DISCRETA_CPP = /usr/bin/g++
DISCRETA_CC = /usr/bin/gcc
CPPFLAGS = -Wall -I$(SRC)/LIB -I$(SRC)/LIB/DISCRETA
LIB = \
	$(SRC)/LIB/TOP_LEVEL/libtoplevel.a \
	$(SRC)/LIB/DISCRETA/libdiscreta.a \
	$(SRC)/LIB/SNAKES_AND_LADDERS/libsnakesandladders.a  \
	$(SRC)/LIB/ACTION/libaction.a \
	$(SRC)/LIB/INCIDENCE/libincidence.a \
	$(SRC)/LIB/GALOIS/libgalois.a 
LFLAGS = -lm



all: test.out test2.out test3.out test4.out \
		partitions.out \
		conic.out hamming.out lehmer.out sarnak.out winnie_li.out \
		field.out design.out change_1_dollar.out \
		puzzle.out \
		pentomino_5x5.out \
		hadamard.out \
		had.out \
		grassmann.out \
		paley.out \
		fibonacci.out \
		johnson.out \
		sequences.out \
		johnson_table.out \
		schlaefli.out \
		all_k_subsets.out \
		cayley_sym_n.out \
		cayley.out \
		shrikhande.out  \
		burnside.out \
		group_ring.out \
		matrix_rank.out \
		alphabet.out \
		alphabet2.out \
		subsets.out \
		orthogonal.out \
		srg.out \
		conjugacy_classes_sym_n.out \
		unrank.out \
		random.out \
		nondecreasing.out
 
#test_group.out 

test.out:  test.o;
	$(DISCRETA_CPP)  test.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

test2.out:  test2.o;
	$(DISCRETA_CPP)  test2.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

test3.out:  test3.o;
	$(DISCRETA_CPP)  test3.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

test4.out:  test4.o;
	$(DISCRETA_CPP)  test4.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

test_group.out:  test_group.o;
	$(DISCRETA_CPP)  test_group.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

partitions.out: partitions.C
	g++ partitions.C -o partitions.out

conic.out:  conic.o;
	$(DISCRETA_CPP)  conic.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

hamming.out:  hamming.o;
	$(DISCRETA_CPP)  hamming.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

lehmer.out:  lehmer.o;
	$(DISCRETA_CC)  lehmer.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

sarnak.out:  sarnak.o;
	$(DISCRETA_CPP)  sarnak.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

winnie_li.out:  winnie_li.o;
	$(DISCRETA_CPP)  winnie_li.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

field.out:  field.o;
	$(DISCRETA_CPP)  field.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

design.out:  design.o;
	$(DISCRETA_CPP)  design.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

change_1_dollar.out:  change_1_dollar.o;
	$(DISCRETA_CPP)  change_1_dollar.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

puzzle.out:  puzzle.o;
	$(DISCRETA_CPP)  puzzle.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

pentomino_5x5.out:  pentomino_5x5.o;
	$(DISCRETA_CPP)  pentomino_5x5.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

hadamard.out:  hadamard.o;
	$(DISCRETA_CPP)  hadamard.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

had.out:  had.o;
	$(DISCRETA_CPP)  had.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

grassmann.out:  grassmann.o;
	$(DISCRETA_CPP)  grassmann.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

paley.out:  paley.o;
	$(DISCRETA_CPP)  paley.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

fibonacci.out:  fibonacci.o;
	$(DISCRETA_CPP)  fibonacci.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

johnson.out:  johnson.o;
	$(DISCRETA_CPP)  johnson.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

sequences.out:  sequences.o;
	$(DISCRETA_CPP)  sequences.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

johnson_table.out:  johnson_table.o;
	$(DISCRETA_CPP)  johnson_table.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

schlaefli.out:  schlaefli.o;
	$(DISCRETA_CPP)  schlaefli.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

all_k_subsets.out:  all_k_subsets.o;
	$(DISCRETA_CPP)  all_k_subsets.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

cayley_sym_n.out:  cayley_sym_n.o;
	$(DISCRETA_CPP)  cayley_sym_n.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

cayley.out:  cayley.o;
	$(DISCRETA_CPP)  cayley.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

shrikhande.out:  shrikhande.o;
	$(DISCRETA_CPP)  shrikhande.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

burnside.out:  burnside.o;
	$(DISCRETA_CPP)  burnside.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

group_ring.out:  group_ring.o;
	$(DISCRETA_CPP)  group_ring.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

matrix_rank.out:  matrix_rank.o;
	$(DISCRETA_CPP)  matrix_rank.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

alphabet.out:  alphabet.o;
	$(DISCRETA_CPP)  alphabet.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

alphabet2.out:  alphabet2.o;
	$(DISCRETA_CPP)  alphabet2.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

subsets.out:  subsets.o;
	$(DISCRETA_CPP)  subsets.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

orthogonal.out:  orthogonal.o;
	$(DISCRETA_CPP)  orthogonal.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

srg.out:  srg.o;
	$(DISCRETA_CPP)  srg.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

conjugacy_classes_sym_n.out:  conjugacy_classes_sym_n.o;
	$(DISCRETA_CPP)  conjugacy_classes_sym_n.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

unrank.out:  unrank.o;
	$(DISCRETA_CPP)  unrank.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

random.out:  random.o;
	$(DISCRETA_CPP)  random.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@

nondecreasing.out:  nondecreasing.o;
	$(DISCRETA_CPP)  nondecreasing.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@



clean:
	- rm *.o *.a *.out *~

test.o: test.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c test.C -o $@
test2.o: test2.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c test2.C -o $@
test4.o: test4.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c test4.C -o $@
test_group.o: test_group.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c test_group.C -o $@
conic.o: conic.C; $(DISCRETA_CPP) $(CPPFLAGS) -c conic.C -o $@
hamming.o: hamming.C; $(DISCRETA_CPP) $(CPPFLAGS) -c hamming.C -o $@
sarnak.o: sarnak.C; $(DISCRETA_CPP) $(CPPFLAGS) -c sarnak.C -o $@
winnie_li.o: winnie_li.C; $(DISCRETA_CPP) $(CPPFLAGS) -c winnie_li.C -o $@
field.o: field.C; $(DISCRETA_CPP) $(CPPFLAGS) -c field.C -o $@
design.o: design.C; $(DISCRETA_CPP) $(CPPFLAGS) -c design.C -o $@
change_1_dollar.o: change_1_dollar.C; $(DISCRETA_CPP) $(CPPFLAGS) -c change_1_dollar.C -o $@
puzzle.o: puzzle.C; $(DISCRETA_CPP) $(CPPFLAGS) -c puzzle.C -o $@
pentomino_5x5.o: pentomino_5x5.C; $(DISCRETA_CPP) $(CPPFLAGS) -c pentomino_5x5.C -o $@
hadamard.o: hadamard.C; $(DISCRETA_CPP) $(CPPFLAGS) -c hadamard.C -o $@
had.o: had.C; $(DISCRETA_CPP) $(CPPFLAGS) -c had.C -o $@
grassmann.o: grassmann.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c grassmann.C -o $@
paley.o: paley.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c paley.C -o $@
fibonacci.o: fibonacci.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c fibonacci.C -o $@
johnson.o: johnson.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c johnson.C -o $@
sequences.o: sequences.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c sequences.C -o $@
johnson_table.o: johnson_table.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c johnson_table.C -o $@
schlaefli.o: schlaefli.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c schlaefli.C -o $@
all_k_subsets.o: all_k_subsets.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c all_k_subsets.C -o $@
cayley_sym_n.o: cayley_sym_n.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c cayley_sym_n.C -o $@
cayley.o: cayley.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c cayley.C -o $@
shrikhande.o: shrikhande.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c shrikhande.C -o $@
burnside.o: burnside.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c burnside.C -o $@
group_ring.o: group_ring.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c group_ring.C -o $@
matrix_rank.o: matrix_rank.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c matrix_rank.C -o $@
alphabet.o: alphabet.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c alphabet.C -o $@
alphabet2.o: alphabet2.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c alphabet2.C -o $@
subsets.o: subsets.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c subsets.C -o $@
orthogonal.o: orthogonal.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c orthogonal.C -o $@
srg.o: srg.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c srg.C -o $@
conjugacy_classes_sym_n.o: conjugacy_classes_sym_n.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c conjugacy_classes_sym_n.C -o $@
unrank.o: unrank.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c unrank.C -o $@
random.o: random.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c random.C -o $@
nondecreasing.o: nondecreasing.C ; $(DISCRETA_CPP) $(CPPFLAGS) -c nondecreasing.C -o $@


tex1:
	~/bin/a2tex.out <puzzle.C >~/TEX.15/ORBITER/CODE/puzzle.tex
	./puzzle.out >puzzle_out.txt
	~/bin/a2tex.out <puzzle_out.txt >~/TEX.15/ORBITER/CODE/puzzle_out.tex
	~/bin/a2tex.out <puzzle25_compressed.txt >~/TEX.15/ORBITER/CODE/puzzle25_compressed.tex



	
