#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
using namespace Eigen;
using namespace std;

class ReactionsSolver {
    public:
        int Num_reactions;
        int Num_components;
        const VectorXd N;
        const VectorXd K;
        VectorXd variablesSet;
        const MatrixXd V;
        MatrixXd V_extended;

        ReactionsSolver(int nr, int nc, const VectorXd &N, const VectorXd &K, const MatrixXd &V)
            : Num_reactions(nr), Num_components(nc), N(N), K(K), V(V)
        {
            variablesSet = VectorXd::Zero(nr+nc);
            V_extended = MatrixXd::Zero(nr+nc,nr+nc);
        }
        void init_extV(){
            V_extended.topLeftCorner(Num_reactions, Num_reactions) = MatrixXd::Zero(Num_reactions, Num_reactions);
            V_extended.topRightCorner(Num_reactions, Num_components) = V.transpose();
            V_extended.bottomLeftCorner(Num_components, Num_reactions) = V;
            V_extended.bottomRightCorner(Num_components, Num_components)= MatrixXd::Zero(Num_components, Num_components);
            V_extended *= -1;
        }

    VectorXd function_F() {
        VectorXd F(variablesSet);
        F.head(Num_reactions) = - K.array().log() - (V.transpose() * variablesSet.tail(Num_components)).array();
        F.tail(Num_components) = - N.array() - (V*variablesSet.head(Num_reactions)).array() + variablesSet.tail(Num_components).array().exp() ;
        return F;
    }

    MatrixXd function_dF(){
        MatrixXd dF = V_extended;
        for ( int i = 0; i < Num_components; i++){
            dF(i + Num_reactions,i + Num_reactions) += exp(variablesSet[Num_reactions + i]);
        }
        return dF;
    }

    VectorXd newtons_method() {
        int l = 0;
        double a = 1.0;
        double e = 1e-8;
        VectorXd variablesSetOld(Num_components + Num_reactions);
        VectorXd A(Num_components +  Num_reactions);
        MatrixXd B(Num_components +  Num_reactions, Num_components +  Num_reactions);
        VectorXd C(Num_components + Num_reactions);
        while(true) {
            l+=1;
            variablesSetOld = variablesSet;
            A = function_F();
            B = function_dF();
            C = B.lu().solve(A);
            variablesSet -= a * C;
            variablesSetOld -= variablesSet;
            if (variablesSetOld.norm() < e){
                    return variablesSet.array().exp();
            }
        }
    }
};

int main(){
    const int Num_reactions = 3;
    const int Num_components = 7;
    MatrixXd V(Num_components, Num_reactions);
    V << -1.0, 0.0 ,0.0,
        0.0, -1.0, 0.0,
        0.0, 0.0, -1.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 1.0,
       -1.0, -1.0, -2.0,
        0.0, 0.0, 1.0;
    VectorXd K(Num_reactions);
    K << pow(10, -14), pow(10, -5.928), pow(10, -8.094);

    VectorXd N(Num_components);
    N << 0.0, 0.0, 0.1, 50.0, 0.0, 0.0, 0.0;

    ReactionsSolver t(Num_reactions, Num_components, N, K, V);

    VectorXd x(Num_components + Num_reactions);
    t.init_extV();
    t.Num_reactions = 3;
    x = t.newtons_method();
    cout << x << endl;
    return 0;
}
