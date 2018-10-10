/**
*
* build with line:
* g++ very_important.cpp hello_world.cpp -o run
*/
#include "very_important.h" // "" checks current working directory
void do_something_wrong()
{
	int foo [5] = { 16, 2, 77, 0, 12071 };
	foo[2] = 42/foo[3];
	printf("%i", foo[2]);
}
int main()
{
	print_message();
	//do_something_wrong(); // to show the debugger
}
