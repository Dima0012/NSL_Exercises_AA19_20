#ifndef __Position__
#define __Position__

class Position {

    private:
        double m_x, m_y, m_z;

    public:

        Position() {
            m_x = 0.;
            m_y = 0.;
            m_z = 0.;
        }

        ~Position(){};

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

        void Reset(){
            m_x = 0.;
            m_y = 0.;
            m_z = 0.;
        }

        void Reset(double x, double y, double z){
            m_x = x;
            m_y = y;
            m_z = z;
        }

};


#endif