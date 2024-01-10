Li_Fitter.exe : histo_Li.cc
	g++ -W -Wall -O2 $^ -o $@ `root-config --cflags --libs`
