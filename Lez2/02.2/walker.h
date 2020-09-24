#ifndef __Walker__
#define __Walker__

class Walker {

    private:
        double m_x, m_y, m_z;

    public:

        Walker() {
            m_x = 0.;
            m_y = 0.;
            m_z = 0.;
        }

        ~Walker(){};

        void Set_X(double a){
            m_x = a;
        }

        void Set_Y(double a){
            m_y = a;
        }

        void Set_Z(double a){
            m_z = a;
        }


        double Get_X(){
            return m_x;
        }

        double Get_Y(){
            return m_y;
        }

        double Get_Z(){
            return m_z;
        }

};


#endif