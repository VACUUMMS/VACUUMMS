#include <iostream>

namespace lib
{
int f();
int g(int a);
}

class Animal
{
    public:
        void live() const;
        virtual void speak() const;
};

class Dog: public Animal
{
    public:
         void speak() const final;
};

class Cat: public Animal
{
    public:
         void speak() const final;
};

