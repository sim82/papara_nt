

//
// test of (possibly bad) idea:
// link this to a utility-source file that contains some 'main-like' code contained in main2.
// FIXME: what happend if more than one source has a main2...
//

int main2( int argc, char *argv[] );

int main( int argc, char *argv[] ) {
    return main2(argc, argv);
}