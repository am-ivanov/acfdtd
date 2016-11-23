#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv) {
	if (argc != 2) {
		cerr << "Usage: " << argv[0] << " columns" << endl;
		return 0;
	}
	int col = atoi(argv[1]);
	string sep;
	if (argc == 3) sep = string(argv[3]);
	else sep = string(" ");
	float val;
	int counter = 0;
	while (1) {
		cin.read(reinterpret_cast<char*>(&val), sizeof(float));
		if (cin.eof()) break;
		cout << val << " ";
		++counter;
		if (counter == col) {
			counter = 0;
			cout << endl;
		}
	}
}
