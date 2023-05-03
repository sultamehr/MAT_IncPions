//File: GetIngredient.h
//Brief: Get a histogram or other TObject from a TFile or TDirectory for one of the cross section extraction steps.
//       Does a little error handling for me automatically.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_GETINGREDIENT_H
#define UTIL_GETINGREDIENT_H

//ROOT includes
#include "TDirectoryFile.h"

//c++ includes
#include <string>
#include <exception>

namespace util
{
  template <class TYPE>
  TYPE* GetIngredient(TDirectoryFile& dir, const std::string& ingredient)
  {
    TObject* obj = dir.Get(ingredient.c_str());
    if(obj == nullptr) throw std::runtime_error("Failed to get " + ingredient + " in " + dir.GetName());
  
    auto typed = dynamic_cast<TYPE*>(obj);
    if(typed == nullptr) throw std::runtime_error(std::string("Found ") + obj->GetName() + ", but it's not the right kind of TObject.");
  
    return typed;
  }
                                                                                                                                         
  template <class TYPE>
  TYPE* GetIngredient(TDirectoryFile& dir, const std::string& ingredient, const std::string& prefix)
  {
    return GetIngredient<TYPE>(dir, prefix + "_" + ingredient);
  }
}

#endif //UTIL_GETINGREDIENT_H
