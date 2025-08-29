/* vacuumms_cpp/FTW.cc */

#include <vacuumms/FTW.hh>

#include <iostream>
#include <typeinfo>

FTWX::FTWX()
{
}

void FTWX::add(FTW* f)
{ 
    std::cout << "Adding object of type: " << typeid(*f).name() << "\n";

    std::cout << "dump says: \n" ;
    f->dump();
    std::cout << "\nend dump.\n" ;
    objects.push_back(f);
}

void FTWX::play()
{
    for (const auto* obj : objects)
        obj->dump();
}

FTW::FTW()
{
}

FTW1::FTW1()
{
}

FTW2::FTW2()
{
}

void FTW::dump() const
{
    std::cout << "dumping an FTW\n\n";
}

void FTW1::dump() const
{
    std::cout << "dumping an FTW1\n\n";
}

void FTW2::dump() const
{
    std::cout << "dumping an FTW2\n\n";
}


