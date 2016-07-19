#include <iostream>
#include <vector>

using namespace std;
void f(int&) {};
void g(const int&) {};

int main(int argc, char* argv[])
{
	int a = 1;
	f(a);
	g(1);
	vector<int> b;
	if (b.empty()) cout << "B is empty!" << endl;
	return 0;
}