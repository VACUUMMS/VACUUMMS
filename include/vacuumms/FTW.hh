/* vacuumms/FTW.hh */

#pragma once

#include <vector>

#include <vacuumms/exports.hh>




class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FTW
{
    public:

        FTW();
        virtual void dump() const;
};


// Trampoline class for Base
class PyFTW : public FTW {
public:
    using FTW::FTW; // Inherit constructors
    void dump() const override {
        PYBIND11_OVERRIDE_PURE(
            void,           // Return type
            FTW,           // Parent class
            dump // Function name
        );
    }
};


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FTWX
{
    public:

        FTWX();
        void add(FTW*);
        void play();

        std::vector<FTW*> objects;
};


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FTW1 : public FTW
{
    public:
        
        FTW1();
        void dump() const override;
};


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FTW2 : public FTW
{
    public:
        
        FTW2();
        void dump() const override;
};



