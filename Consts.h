#ifndef REACTIONS_CONSTS_H_
#define REACTIONS_CONSTS_H_

namespace Reactions
{
  class Consts
  {
  public:

    /*A09*/
    static constexpr double k2on          = 1.0e7;
    static constexpr double N2b           = 1000.0;
    static constexpr double N2se          = 1000.0;
    static constexpr double k2off         = 5.9;
    static constexpr double kz2promplus   = 1.03e8;
    static constexpr double kz2promminus  = 1.0;
    /*A10*/
    static constexpr double kz2promcat    = 30.0;
    static constexpr double kz5e2mcat     = 0.23;
    static constexpr double kz5e2mminus   = 1.0;
    static constexpr double kz5e2mplus    = 1.73e7;
    static constexpr double kz8e2mcat     = 0.9;
    static constexpr double kz8e2mminus   = 1.0;
    static constexpr double kz8e2mplus    = 2.64e7;
    /*A11*/
    static constexpr double k5on          = 5.7e7;
    static constexpr double N5b           = 3000;
    static constexpr double N5se          = 3000;
    static constexpr double k5off         = 0.17;
    static constexpr double kz5e10mplus   = 1.0e8;
    static constexpr double kz5e10mminus  = 1.0;
    /*A12*/
    static constexpr double kz5e10mcat    = 4.6e-2;
    static constexpr double kprominus     =  0.01;
    static constexpr double kproplus      = 1.0e8; 
    static constexpr double kapce5mplus   = 1.2e8;
    static constexpr double kapce5mminus   =  1.0;
    /*A13*/
    static constexpr double k8on          = 5.0e7;
    static constexpr double N8b           = 450.0;
    static constexpr double N8se          = 450.0; 
    static constexpr double k8off         = 0.17; 
    static constexpr double kz8e10mplus   = 5.1e7; 
    static constexpr double kz8e10mminus  = 1.0;
    /*A14*/
    static constexpr double kz8e10mcat    = 2.3e-2;
    static constexpr double ktenminus     = 0.01; 
    static constexpr double ktenplus      = 1.0e8; 
    static constexpr double kapce8mplus   = 1.2e8; 
    static constexpr double kapce8mminus  = 1.0;
    /*A15*/
    static constexpr double k9on          = 1.0e7; 
    static constexpr double N9b           = 250.0; 
    static constexpr double N9se          = 250.0;
    static constexpr double k9off         = 2.5e-2;
    /*A17*/
    static constexpr double k10on         = 1.0e7;
    static constexpr double N10b          = 2700.0; 
    static constexpr double N10se         = 2700.0; 
    static constexpr double k10off        = 2.5e-2; 
    static constexpr double kz10tenmminus  = 1.0;
    static constexpr double kz10tenmplus  = 1.31e8;
    /*A18*/
    static constexpr double kz10tenmcat   = 20.0;
    /*A19*/
    /*A27*/
    static constexpr double kapce5mcat    = 0.5;
    /*A28*/
    static constexpr double kapce8mcat    = 0.5;
    /*A29*/
    static constexpr double N9starb       = 250.0;
    static constexpr double N9starse      = 250.0;
  };
}


#endif //REACTIONS_CONSTEXPRS_H_
