#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define FIXOUTVEC(a,b,c,d,e) FIXOUT(a,b << "[   ");a.width(5);a << std::left << c << " , ";a.width(5);a << std::left << d << " , ";a.width(5);a << std::left << e << "   ]" << std::endl;
#define SUCCESS(a) a << "Done" << std::endl;
