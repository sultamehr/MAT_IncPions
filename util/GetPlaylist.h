//File: GetPlaylist.h
//Brief: A function that looks at an AnaTuple file and figures out the playlist name for the first event.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_GETPLAYLIST_H
#define UTIL_GETPLAYLIST_H

//c++ includes
#include <string>

namespace PlotUtils
{
  class TreeWrapper;
}

namespace util
{
  std::string GetPlaylist(PlotUtils::TreeWrapper& tree, const bool isMC);
}

#endif //UTIL_GETPLAYLIST_H
