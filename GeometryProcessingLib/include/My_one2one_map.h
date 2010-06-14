#pragma once

#include <GL/glew.h>
#include <GL/glut.h>

#include <map>

#include "Sort.h"
using namespace std;

namespace Another
{

    /*****************************************************************************************************/
    /*****************This one to one map is used for convert m file to HDS structure************/
    /*****************************************************************************************************/
    
    template<class Key, class Value, class Sort>
    class My_One2OneMap
     {
     
        private:
            map<Key,Value> vertexMapO;
            map<Value, Key,Sort> vertexMapOI;

        public:
            typedef typename map<Key,Value>::iterator MAPIT;
            typedef typename map<Value, Key,Sort>::iterator MAPIIT;	
            
            My_One2OneMap(){}
            ~My_One2OneMap(){}
            
            bool InsertvertexHandleO(Key key, Value vo)
            {
                this->vertexMapO.insert( pair<Key, Value>( key, vo) );
                this->vertexMapOI.insert( pair<Value, Key>( vo, key) );
                return true;
            }

            bool DeletevertexHandleO(Key key)
            {
                MAPIT it;
                it = vertexMapO.find(key);
                vertexMapO.erase( it );
                return true;        
            }	

            Key FindIndexOI(Value vo)
            {
                MAPIIT it = vertexMapOI.find(vo);
                if( it == vertexMapOI.end() )  
                {
                    return -1;
                }
                else
                 {
                    return it->second;
                 }
             }

            Value FindVertexhandleO(Key key)
             {
                MAPIT it = vertexMapO.find(key);
                if( it == vertexMapO.end() )
                {
                    return NULL;
                }
                else
                {
                    return it->second;
                }
            }
            
            MAPIT begin()
            {
                return vertexMapO.begin();
            }
            
            MAPIT end()
            {
                return vertexMapO.end();
            }
        };

    }