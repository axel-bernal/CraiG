///////////////////////////////////////////////////////////////////////////////
//
// ObjectFactory
//
// The ObjectFactory class is a object factory implementation.  It allows users
// to register and unregister classes during run-time by specifying a
// user-defined unique identifier per class.  Instances of any registered class
// can then be instantiated simply by calling the create method and passing the
// proper unique identifier.
//
////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2004 Robert Geiman.
//
// Permission to copy, modify, and use this code for personal and commercial
// software is granted provided this copyright notice appears in all copies.
// This software is provided "as is" without any expressed or implied warranty.
//
// Any comments or questions can be sent to: rgeiman@buckeye-express.com
//
///////////////////////////////////////////////////////////////////////////////


#ifndef OBJECT_FACTORY_H
#define OBJECT_FACTORY_H

#include <map>


namespace lless {

    template<typename Type>
        struct Type2Type
        {
            typedef Type OriginalType;
        };


    /////////////////////////////////////////////////////////////////////////////
    // Because some compilers, such as Visual C++ 6.0, do not support partial
    // template specialization we must create seperate ObjectFactory classes, one
    // for each constructor parameter we allow the user to use.
    /////////////////////////////////////////////////////////////////////////////

    /**
     * In general, an ObjectFactoryN class is a object factory implementation.
     * It allows users to register and unregister classes with N-parameter
     * constructors during run-time by specifying a user-defined unique
     * identifier per class.  Instances of any registered class can then be
     * instantiated simply by calling the create method and passing the proper
     * unique identifier.
     ************************************************************************/

    template<typename BaseClassType, typename Param1Type, typename Param2Type, typename UniqueIdType> class ObjectFactory2 {
    protected:
        class CreateObjectBase
        {
        public:
            virtual ~CreateObjectBase() {}
            virtual BaseClassType *operator()(Param1Type a1, Param2Type a2) = 0;
        };
        template<typename ClassType> class CreateObject : public CreateObjectBase
        {
        public:
            BaseClassType *operator()(Param1Type a1, Param2Type a2)
            {
                return (BaseClassType *)(new ClassType(a1, a2));
            }
        };
    public:
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::const_iterator ConstIterator;
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::iterator Iterator;
        ~ObjectFactory2() {
            Iterator iter = m_object_creator.begin();
            while (iter != m_object_creator.end()) {
                delete (*iter).second;
                ++iter;
            }
        }

        template<typename ClassType> bool Register(UniqueIdType unique_id, Type2Type<ClassType>) {
            if (m_object_creator.find(unique_id) != m_object_creator.end())
                return false;
            m_object_creator[unique_id] = new CreateObject<ClassType>;
            return true;
        }

        bool Unregister(UniqueIdType unique_id) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter != m_object_creator.end()) {
                delete (*iter).second;
                m_object_creator.erase(iter);
                return true;
            }
            return false;
        }

        BaseClassType *Create(UniqueIdType unique_id, Param1Type a1, Param2Type a2) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter == m_object_creator.end())
                return NULL;
            return ((*iter).second)->operator()(a1, a2);
        }

        ConstIterator GetBegin() const {
            return m_object_creator.begin();
        }

        Iterator GetBegin() {
            return m_object_creator.begin();
        }

        ConstIterator GetEnd() const {
            return m_object_creator.end();
        }

        Iterator GetEnd() {
            return m_object_creator.end();
        }

    protected:
        std::map<UniqueIdType, CreateObjectBase*> m_object_creator;
    };



    /*
     * The ObjectFactory3 class is a object factory implementation.
     * It allows users to register and unregister classes with 3-parameter
     * constructors during run-time by specifying a user-defined unique
     * identifier per class.  Instances of any registered class can then be
     * instantiated simply by calling the create method and passing the proper
     * unique identifier.
     ************************************************************************/

    template<typename BaseClassType, typename Param1Type, typename Param2Type, typename Param3Type, typename UniqueIdType> class ObjectFactory3 {
    protected:
        class CreateObjectBase
        {
        public:
            virtual ~CreateObjectBase() {}
            virtual BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3) = 0;
        };
        template<typename ClassType> class CreateObject : public CreateObjectBase
        {
        public:
            BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3)
            {
                return (BaseClassType *)(new ClassType(a1, a2, a3));
            }
        };
    public:
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::const_iterator ConstIterator;
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::iterator Iterator;
        ~ObjectFactory3() {
            Iterator iter = m_object_creator.begin();
            while (iter != m_object_creator.end()) {
                delete (*iter).second;
                ++iter;
            }
        }

        template<typename ClassType> bool Register(UniqueIdType unique_id, Type2Type<ClassType>) {
            if (m_object_creator.find(unique_id) != m_object_creator.end())
                return false;
            m_object_creator[unique_id] = new CreateObject<ClassType>;
            return true;
        }

        bool Unregister(UniqueIdType unique_id) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter != m_object_creator.end()) {
                delete (*iter).second;
                m_object_creator.erase(iter);
                return true;
            }
            return false;
        }

        BaseClassType *Create(UniqueIdType unique_id, Param1Type a1, Param2Type a2, Param3Type a3) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter == m_object_creator.end())
                return NULL;
            return ((*iter).second)->operator()(a1, a2, a3);
        }

        ConstIterator GetBegin() const {
            return m_object_creator.begin();
        }

        Iterator GetBegin() {
            return m_object_creator.begin();
        }

        ConstIterator GetEnd() const {
            return m_object_creator.end();
        }

        Iterator GetEnd() {
            return m_object_creator.end();
        }

    protected:
        std::map<UniqueIdType, CreateObjectBase*> m_object_creator;
    };



    /*
     * The ObjectFactory4 class is a object factory implementation.
     * It allows users to register and unregister classes with 4-parameter
     * constructors during run-time by specifying a user-defined unique
     * identifier per class.  Instances of any registered class can then be
     * instantiated simply by calling the create method and passing the proper
     * unique identifier.
     ************************************************************************/
    template<typename BaseClassType, typename Param1Type, typename Param2Type, typename Param3Type, typename Param4Type, typename UniqueIdType> class ObjectFactory4 {
    protected:
        class CreateObjectBase
        {
        public:
            virtual ~CreateObjectBase() {}
            virtual BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4) = 0;
        };
        template<typename ClassType> class CreateObject : public CreateObjectBase
        {
        public:
            BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4)
            {
                return (BaseClassType *)(new ClassType(a1, a2, a3, a4));
            }
        };
    public:
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::const_iterator ConstIterator;
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::iterator Iterator;
        ~ObjectFactory4() {
            Iterator iter = m_object_creator.begin();
            while (iter != m_object_creator.end()) {
                delete (*iter).second;
                ++iter;
            }
        }

        template<typename ClassType> bool Register(UniqueIdType unique_id, Type2Type<ClassType>) {
            if (m_object_creator.find(unique_id) != m_object_creator.end())
                return false;
            m_object_creator[unique_id] = new CreateObject<ClassType>;
            return true;
        }

        bool Unregister(UniqueIdType unique_id) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter != m_object_creator.end()) {
                delete (*iter).second;
                m_object_creator.erase(iter);
                return true;
            }
            return false;
        }

        BaseClassType *Create(UniqueIdType unique_id, Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter == m_object_creator.end())
                return NULL;
            return ((*iter).second)->operator()(a1, a2, a3, a4);
        }

        ConstIterator GetBegin() const {
            return m_object_creator.begin();
        }

        Iterator GetBegin() {
            return m_object_creator.begin();
        }

        ConstIterator GetEnd() const {
            return m_object_creator.end();
        }

        Iterator GetEnd() {
            return m_object_creator.end();
        }

    protected:
        std::map<UniqueIdType, CreateObjectBase*> m_object_creator;
    };



    /*
     * The ObjectFactory5 class is a object factory implementation.
     * It allows users to register and unregister classes with 5-parameter
     * constructors during run-time by specifying a user-defined unique
     * identifier per class.  Instances of any registered class can then be
     * instantiated simply by calling the create method and passing the proper
     * unique identifier.
     ************************************************************************/

    template<typename BaseClassType, typename Param1Type, typename Param2Type, typename Param3Type, typename Param4Type, typename Param5Type, typename UniqueIdType> class ObjectFactory5 {
    protected:
        class CreateObjectBase
        {
        public:
            virtual ~CreateObjectBase() {}
            virtual BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5) = 0;
        };
        template<typename ClassType> class CreateObject : public CreateObjectBase
        {
        public:
            BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5)
            {
                return (BaseClassType *)(new ClassType(a1, a2, a3, a4, a5));
            }
        };
    public:
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::const_iterator ConstIterator;
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::iterator Iterator;
        ~ObjectFactory5() {
            Iterator iter = m_object_creator.begin();
            while (iter != m_object_creator.end()) {
                delete (*iter).second;
                ++iter;
            }
        }

        template<typename ClassType> bool Register(UniqueIdType unique_id, Type2Type<ClassType>) {
            if (m_object_creator.find(unique_id) != m_object_creator.end())
                return false;
            m_object_creator[unique_id] = new CreateObject<ClassType>;
            return true;
        }

        bool Unregister(UniqueIdType unique_id) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter != m_object_creator.end()) {
                delete (*iter).second;
                m_object_creator.erase(iter);
                return true;
            }
            return false;
        }

        BaseClassType *Create(UniqueIdType unique_id, Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter == m_object_creator.end())
                return NULL;
            return ((*iter).second)->operator()(a1, a2, a3, a4, a5);
        }

        ConstIterator GetBegin() const {
            return m_object_creator.begin();
        }

        Iterator GetBegin() {
            return m_object_creator.begin();
        }

        ConstIterator GetEnd() const {
            return m_object_creator.end();
        }

        Iterator GetEnd() {
            return m_object_creator.end();
        }

    protected:
        std::map<UniqueIdType, CreateObjectBase*> m_object_creator;
    };



    /*
     * The ObjectFactory6 class is a object factory implementation.
     * It allows users to register and unregister classes with 6-parameter
     * constructors during run-time by specifying a user-defined unique
     * identifier per class.  Instances of any registered class can then be
     * instantiated simply by calling the create method and passing the proper
     * unique identifier.
     ************************************************************************/

    template<typename BaseClassType, typename Param1Type, typename Param2Type, typename Param3Type, typename Param4Type, typename Param5Type, typename Param6Type, typename UniqueIdType> class ObjectFactory6 {
    protected:
        class CreateObjectBase
        {
        public:
            virtual ~CreateObjectBase() {}
            virtual BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5, Param6Type a6) = 0;
        };
        template<typename ClassType> class CreateObject : public CreateObjectBase
        {
        public:
            BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5, Param6Type a6)
            {
                return (BaseClassType *)(new ClassType(a1, a2, a3, a4, a5, a6));
            }
        };
    public:
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::const_iterator ConstIterator;
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::iterator Iterator;
        ~ObjectFactory6() {
            Iterator iter = m_object_creator.begin();
            while (iter != m_object_creator.end()) {
                delete (*iter).second;
                ++iter;
            }
        }

        template<typename ClassType> bool Register(UniqueIdType unique_id, Type2Type<ClassType>) {
            if (m_object_creator.find(unique_id) != m_object_creator.end())
                return false;
            m_object_creator[unique_id] = new CreateObject<ClassType>;
            return true;
        }

        bool Unregister(UniqueIdType unique_id) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter != m_object_creator.end()) {
                delete (*iter).second;
                m_object_creator.erase(iter);
                return true;
            }
            return false;
        }

        BaseClassType *Create(UniqueIdType unique_id, Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5, Param6Type a6) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter == m_object_creator.end())
                return NULL;
            return ((*iter).second)->operator()(a1, a2, a3, a4, a5, a6);
        }

        ConstIterator GetBegin() const {
            return m_object_creator.begin();
        }

        Iterator GetBegin() {
            return m_object_creator.begin();
        }

        ConstIterator GetEnd() const {
            return m_object_creator.end();
        }

        Iterator GetEnd() {
            return m_object_creator.end();
        }

    protected:
        std::map<UniqueIdType, CreateObjectBase*> m_object_creator;
    };



    /*
     * The ObjectFactory7 class is a object factory implementation.
     * It allows users to register and unregister classes with 7-parameter
     * constructors during run-time by specifying a user-defined unique
     * identifier per class.  Instances of any registered class can then be
     * instantiated simply by calling the create method and passing the proper
     * unique identifier.
     ************************************************************************/

    template<typename BaseClassType, typename Param1Type, typename Param2Type, typename Param3Type, typename Param4Type, typename Param5Type, typename Param6Type, typename Param7Type, typename UniqueIdType> class ObjectFactory7 {
    protected:
        class CreateObjectBase
        {
        public:
            virtual ~CreateObjectBase() {}
            virtual BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5, Param6Type a6, Param7Type a7) = 0;
        };
        template<typename ClassType> class CreateObject : public CreateObjectBase
        {
        public:
            BaseClassType *operator()(Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5, Param6Type a6, Param7Type a7)
            {
                return (BaseClassType *)(new ClassType(a1, a2, a3, a4, a5, a6, a7));
            }
        };
    public:
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::const_iterator ConstIterator;
        typedef typename std::map<UniqueIdType, CreateObjectBase*>::iterator Iterator;
        ~ObjectFactory7() {
            Iterator iter = m_object_creator.begin();
            while (iter != m_object_creator.end()) {
                delete (*iter).second;
                ++iter;
            }
        }

        template<typename ClassType> bool Register(UniqueIdType unique_id, Type2Type<ClassType>) {
            if (m_object_creator.find(unique_id) != m_object_creator.end())
                return false;
            m_object_creator[unique_id] = new CreateObject<ClassType>;
            return true;
        }

        bool Unregister(UniqueIdType unique_id) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter != m_object_creator.end()) {
                delete (*iter).second;
                m_object_creator.erase(iter);
                return true;
            }
            return false;
        }

        BaseClassType *Create(UniqueIdType unique_id, Param1Type a1, Param2Type a2, Param3Type a3, Param4Type a4, Param5Type a5, Param6Type a6, Param7Type a7) {
            Iterator iter = m_object_creator.find(unique_id);
            if (iter == m_object_creator.end())
                return NULL;
            return ((*iter).second)->operator()(a1, a2, a3, a4, a5, a6, a7);
        }

        ConstIterator GetBegin() const {
            return m_object_creator.begin();
        }

        Iterator GetBegin() {
            return m_object_creator.begin();
        }

        ConstIterator GetEnd() const {
            return m_object_creator.end();
        }

        Iterator GetEnd() {
            return m_object_creator.end();
        }

    protected:
        std::map<UniqueIdType, CreateObjectBase*> m_object_creator;
    };

}

#endif
