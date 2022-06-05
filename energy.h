#pragma once

class Energy{
public:
    PS::F64 etot;
    PS::F64 ekin;
    PS::F64 ephi_sun;
    PS::F64 ephi_planet;
    PS::F64 ephi;
    PS::F64 ephi_d;
    PS::F64 edisp;

    

    Energy(){ etot = ekin = ephi_sun = ephi_planet = ephi = ephi_d = edisp = 0.; }
    
    Energy(PS::F64 ekin0,
           PS::F64 ephi_sun0,
           PS::F64 ephi0,
           PS::F64 ephi_d0,
           PS::F64 edisp0){
        ekin = ekin0;
        ephi_sun = ephi_sun0;
        ephi = ephi0;
        ephi_d = ephi_d0;
        edisp = edisp0;
        ephi_planet = ephi + ephi_d;
        etot = ekin + ephi_sun + ephi_planet;
    }

    template<class Tpsys>
    void calcEnergy(const Tpsys & pp,
                    const bool clear=true){
        if ( clear ) etot = ekin = ephi_sun = ephi_planet = ephi = ephi_d = 0.0;

        PS::F64 ekin_loc = 0.0;
        PS::F64 ephi_sun_loc = 0.0;
        PS::F64 ephi_loc = 0.0;
        PS::F64 ephi_d_loc = 0.0;
        
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < n_loc; i++){
            ekin_loc     += pp[i].mass * pp[i].vel * pp[i].vel;
            ephi_sun_loc += pp[i].mass * pp[i].phi_s;
#ifndef CORRECT_NEIGHBOR
            ephi_loc     += pp[i].mass * pp[i].phi;
#else
            ephi_loc     += pp[i].mass * (pp[i].phi + pp[i].phi_correct);
#endif
            ephi_d_loc   += pp[i].mass * pp[i].phi_d;
        }
        ekin_loc *= 0.5;
        ephi_loc *= 0.5;
        ephi_d_loc *= 0.5;

        ekin     += PS::Comm::getSum(ekin_loc);
#ifdef INDIRECT_TERM
        ekin     += getIndirectEnergy(pp);
#endif
        ephi_sun += PS::Comm::getSum(ephi_sun_loc);
        ephi     += PS::Comm::getSum(ephi_loc);
        ephi_d   += PS::Comm::getSum(ephi_d_loc);
        ephi_planet =  ephi + ephi_d;
        etot = ekin + ephi_sun + ephi_planet;
    }

    template<class Tpsys>
    void calcEnergy_output(const Tpsys & pp,
                           PS::F64 time_sys,
                           std::ofstream & fp,
                           bool flag){
        bool clear = true;
        if ( clear ) etot = ekin = ephi_sun = ephi_planet = ephi = ephi_d = 0.0;

        PS::F64 ekin_loc = 0.0;
        PS::F64 ephi_sun_loc = 0.0;
        PS::F64 ephi_loc = 0.0;
        PS::F64 ephi_d_loc = 0.0;
        
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < n_loc; i++){
            ekin_loc     += pp[i].mass * pp[i].vel * pp[i].vel;
            ephi_sun_loc += pp[i].mass * pp[i].phi_s;
#ifndef CORRECT_NEIGHBOR
            ephi_loc     += pp[i].mass * pp[i].phi;
#else
            ephi_loc     += pp[i].mass * (pp[i].phi + pp[i].phi_correct);
#endif
            ephi_d_loc   += pp[i].mass * pp[i].phi_d;
        }
        ekin_loc *= 0.5;
        ephi_loc *= 0.5;
        ephi_d_loc *= 0.5;

        ekin     += PS::Comm::getSum(ekin_loc);
#ifdef INDIRECT_TERM
        ekin     += getIndirectEnergy(pp);
#endif
        ephi_sun += PS::Comm::getSum(ephi_sun_loc);
        ephi     += PS::Comm::getSum(ephi_loc);
        ephi_d   += PS::Comm::getSum(ephi_d_loc);
        ephi_planet =  ephi + ephi_d;
        etot = ekin + ephi_sun + ephi_planet;
        if(flag)
        {
            fp << std::scientific<<std::setprecision(8)<<"after crossing OMF"<< "\t" <<time_sys << "\t"
               << ekin << "\t" << ephi_sun << "\t" << ephi  << "\t"
               << ephi_d << "\t" << ephi_planet << "\t"
               << etot<< std::endl;
        }
        else
        {
            fp << std::scientific<<std::setprecision(8)<<"before crossing OMF"<< "\t" <<time_sys << "\t"
               << ekin << "\t" << ephi_sun << "\t" << ephi  << "\t"
               << ephi_d << "\t" << ephi_planet << "\t"
               << etot<< std::endl;
        }
        for(PS::S32 i = 0; i < n_loc; i++)
        {
            fp  << std::scientific<<std::setprecision(16)
                << time_sys << "\t" << pp[i].id << "\t" << pp[i].mass  << "\t"
                << pp[i].pos.x << "\t" << pp[i].pos.y  << "\t" << pp[i].pos.z <<"\t"
                << pp[i].vel.x << "\t" << pp[i].vel.y  << "\t" << pp[i].vel.z <<"\t"
                <<0.5*pp[i].mass*(pp[i].vel.x*pp[i].vel.x+pp[i].vel.y*pp[i].vel.y)<<"\t"
                <<sqrt(pp[i].pos.x*pp[i].pos.x+pp[i].pos.y*pp[i].pos.y)<<"\t"
                << pp[i].phi << "\t" << pp[i].phi_d << "\t" << pp[i].phi_s<< "\t"
                << std::endl;
        }
    }
    
    template<class Tpsys>
    void calcEnergy_output(const Tpsys & pp,
                           PS::F64 time_sys,
                           std::ofstream & fp){
        bool clear = true;
        if ( clear )
        {
            etot = ekin = ephi_sun = ephi_planet = ephi = ephi_d = 0.0;
        } 

        PS::F64 ekin_loc = 0.0;
        PS::F64 ephi_sun_loc = 0.0;
        PS::F64 ephi_loc = 0.0;
        PS::F64 ephi_d_loc = 0.0;
        PS::F64 etot_ideal = 0.0;
        PS::F64 ekin_ideal = 0.0;
        PS::F64 ephi_sun_ideal = 0.0;
        PS::F64 ephi_ideal = 0.0;
        const PS::F64 m_sun = FP_t::m_sun;
        const PS::F64 eps2  = FP_t::eps2_sun;
        const PS::S32 n_loc = pp.getNumberOfParticleLocal();

        for(PS::S32 i = 0; i < n_loc; i++){
            ekin_loc     += pp[i].mass * pp[i].vel * pp[i].vel;
            ephi_sun_loc += pp[i].mass * pp[i].phi_s;
#ifndef CORRECT_NEIGHBOR
            ephi_loc     += pp[i].mass * pp[i].phi;
#else
            ephi_loc     += pp[i].mass * (pp[i].phi + pp[i].phi_correct);
#endif
            ephi_d_loc   += pp[i].mass * pp[i].phi_d;
        }
        ekin_loc *= 0.5;
        ephi_loc *= 0.5;
        ephi_d_loc *= 0.5;

        ekin     += PS::Comm::getSum(ekin_loc);
#ifdef INDIRECT_TERM
        ekin     += getIndirectEnergy(pp);
#endif
        ephi_sun += PS::Comm::getSum(ephi_sun_loc);
        ephi     += PS::Comm::getSum(ephi_loc);
        ephi_d   += PS::Comm::getSum(ephi_d_loc);
        ephi_planet =  ephi + ephi_d;
        etot = ekin + ephi_sun + ephi_planet;
/***************:ここから位置と速度自体から出した全エネルギーを計算****************/
        for(PS::S32 i = 0; i < n_loc; i++){
            PS::F64vec posi = pp[i].pos;
            PS::F64vec dr_i = - posi; 
            PS::F64    r2inv_i = 1. / (dr_i*dr_i + eps2);
            PS::F64    rinv_i  = sqrt(r2inv_i);
            ekin_ideal     += 0.5*pp[i].mass * pp[i].vel * pp[i].vel;
            ephi_sun_ideal -= pp[i].mass * m_sun * rinv_i;
            for(PS::S32 j = i+1; j < n_loc; j++)
            {
                PS::F64vec r_ij = pp[j].pos - pp[i].pos; 
                PS::F64 dr_ij = sqrt(r_ij*r_ij);
                PS::F64 r2inv_ij = 1. / (dr_ij*dr_ij+eps2);
                PS::F64 rinv_ij = sqrt(r2inv_ij); 
#ifndef CORRECT_NEIGHBOR
                ephi_ideal     -= pp[i].mass * pp[j].mass * rinv_ij;
#else
                ephi_idael     += pp[i].mass * (pp[i].phi + pp[i].phi_correct);
#endif
            }
        }
        etot_ideal = ekin_ideal + ephi_sun_ideal + ephi_ideal;
        fp  <<std::fixed << time_sys << "\t"<<std::scientific<<std::setprecision(8)
            << ekin << "\t" << ephi_sun << "\t" << ephi  << "\t"
            << ephi_d << "\t" << ephi_planet << "\t"
            << etot<< "\t" <<std::endl;
        fp  <<std::fixed << time_sys << "\t"<<std::scientific<<std::setprecision(8)
            << ekin_ideal << "\t" << ephi_sun_ideal << "\t" << ephi_ideal  << "\t"
            << etot_ideal<< "\t" <<std::endl;
    }
    
    PS::F64 calcEnergyError(const Energy e_init){
        //return (etot - e_init.etot - edisp)/e_init.etot;
        return (etot - e_init.etot - edisp)/etot;   //散逸エネルギーのものを引くとかなり小さな値になるはずで、そうなるとエネルギー保存がしっかりされていることがわかる(大体10^-8くらいであればよい)
    }

    PS::F64 calcEnergyError(const Energy e_init,
                            PS::F64 time_sys,
                            std::ofstream & fp,
                            PS::F64 dekin_d){
        //return (etot - e_init.etot - edisp)/e_init.etot;
        fp  <<std::fixed << time_sys << "\t"<<std::scientific<<std::setprecision(16)
            << ekin << "\t" << ephi_sun << "\t" << ephi  << "\t"
            << ephi_d << "\t" << e_init.etot << "\t" << edisp << "\t" <<dekin_d <<"\t"
            << etot << "\t" << (etot - edisp)
            << std::endl;
        return (etot - e_init.etot - edisp)/etot;   //散逸エネルギーのものを引くとかなり小さな値になるはずで、そうなるとエネルギー保存がしっかりされていることがわかる(大体10^-8くらいであればよい)
    }

    template<class Tpsys>
    void check_gas_drag_energy_change(const Tpsys & pp,
                                      PS::F64 dekin_d,
                                      PS::F64 edisp_gd,
                                      PS::F64 time_sys,
                                      std::ofstream & fp,
                                      bool flag){
            const PS::S32 n_loc = pp.getNumberOfParticleLocal();
            /*fp  <<std::fixed<< time_sys << "\t"<<std::scientific<<std::setprecision(16)
                << dekin_d << "\t" << edisp_gd
                << std::endl;*/
            PS::F64 ekin_tot = 0.0;
            if(flag)
            {
                fp  <<"Before vel kick"<<"\t"<<std::fixed<< time_sys << "\t"<<std::scientific<<std::setprecision(16)
                << "\t" << edisp_gd << "\t" << edisp << std::endl;
                for(PS::S32 i = 0; i < n_loc; i++)
                {
                    fp  << std::scientific<<std::setprecision(12)
                    << pp[i].id << "\t" << pp[i].mass  << "\t"
                    << pp[i].pos.x << "\t" << pp[i].pos.y  << "\t"
                    << pp[i].vel.x << "\t" << pp[i].vel.y  << "\t"
                    << pp[i].acc.x << "\t" << pp[i].acc.y  << "\t"
                    << pp[i].acc_gd.x << "\t" << pp[i].acc_gd.y<<"\t"
                    <<0.5*pp[i].mass*(pp[i].vel.x*pp[i].vel.x+pp[i].vel.y*pp[i].vel.y)<<"\t"
                    <<sqrt(pp[i].pos.x*pp[i].pos.x+pp[i].pos.y*pp[i].pos.y)<<"\t"<<std::endl;
                    ekin_tot += 0.5*pp[i].mass*(pp[i].vel.x*pp[i].vel.x+pp[i].vel.y*pp[i].vel.y);
                }   
                fp<<ekin_tot<<std::endl;
            }
            else
            {   
                fp  << "After vel kick" << "\t"<<std::scientific<<std::setprecision(16)
                << "\t" << dekin_d << "\t" << edisp_gd << "\t" << edisp
                << std::endl;
                for(PS::S32 i = 0; i < n_loc; i++)
                {
                    fp  << std::scientific<<std::setprecision(12)
                    << pp[i].id << "\t" << pp[i].mass  << "\t"
                    << pp[i].pos.x << "\t" << pp[i].pos.y  << "\t"
                    << pp[i].vel.x << "\t" << pp[i].vel.y  << "\t"
                    << pp[i].acc.x << "\t" << pp[i].acc.y  << "\t"
                    << pp[i].acc_gd.x << "\t" << pp[i].acc_gd.y  << "\t"
                    <<0.5*pp[i].mass*(pp[i].vel.x*pp[i].vel.x+pp[i].vel.y*pp[i].vel.y)<<"\t"
                    <<sqrt(pp[i].pos.x*pp[i].pos.x+pp[i].pos.y*pp[i].pos.y)<<"\t"<<std::endl;
                    ekin_tot += 0.5*pp[i].mass*(pp[i].vel.x*pp[i].vel.x+pp[i].vel.y*pp[i].vel.y);
                }   
                fp<<ekin_tot<<std::endl;
            }
             
        }
        template<class Tpsys>
        void check_gas_drag_energy_change(const Tpsys & pp,
                                      PS::F64 time_sys,
                                      std::ofstream & fp,
                                      bool flag){
            const PS::S32 n_loc = pp.getNumberOfParticleLocal();
            fp  <<std::fixed<< time_sys
                << std::endl;
            if(flag)
            {
                for(PS::S32 i = 0; i < n_loc; i++)
                {
                    fp  << std::scientific<<std::setprecision(12)
                    << pp[i].id << "\t" << pp[i].mass  << "\t"
                    << pp[i].pos.x << "\t" << pp[i].pos.y  << "\t" << pp[i].pos.z <<"\t"
                    << pp[i].vel.x << "\t" << pp[i].vel.y  << "\t" << pp[i].vel.z <<"\t"
                    << pp[i].acc_gd.x << "\t" << pp[i].acc_gd.y  << "\t" << pp[i].acc_gd.z <<"\t"
                    <<0.5*pp[i].mass*(pp[i].vel.x*pp[i].vel.x+pp[i].vel.y*pp[i].vel.y)<<"\t"
                    <<sqrt(pp[i].pos.x*pp[i].pos.x+pp[i].pos.y*pp[i].pos.y)
                    << std::endl;
                }   
            } 
        }
    };


class FileHeader{
public:
    PS::S32 n_body;
    PS::S32 id_next;
    PS::F64 time;
    //PS::F64 etot0;
    //PS::F64 etot1;
    //PS::F64 edisp;
    Energy e_init;
    Energy e_now;

    void copy(const FileHeader fh) {
        n_body  = fh.n_body;
        id_next = fh.id_next;
        time    = fh.time;
        e_init  = fh.e_init;
        e_now   = fh.e_now;
    }
    
    PS::S32 readAscii(FILE * fp) {
        if ( !fscanf(fp, "%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                     &time, &n_body, &id_next,
                     &e_init.etot, &e_init.ekin, &e_init.ephi_sun, &e_init.ephi_planet, &e_init.edisp,
                     &e_now.etot,  &e_now.ekin,  &e_now.ephi_sun,  &e_now.ephi_planet,  &e_now.edisp) ) {
            
            errorMessage("The header has NOT been correctly read.");
            PS::Abort();
        }
        return n_body;
    }
    void writeAscii(FILE* fp) const {
        if ( !fprintf(fp,
                      "%g\t%d\t%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\n",
                      time, n_body, id_next,
                      e_init.etot, e_init.ekin, e_init.ephi_sun, e_init.ephi_planet, e_init.edisp,
                      e_now.etot,  e_now.ekin,  e_now.ephi_sun,  e_now.ephi_planet,  e_now.edisp) ) {
            
            errorMessage("The header has NOT been correctly written.");
            PS::Abort();
        }
    }
    PS::S32 readBinary(FILE * fp) {
        FileHeader buf;
        if ( !fread(&buf, sizeof(buf), 1, fp) ) {
            errorMessage("The header has NOT been correctly read.");
            PS::Abort();
        }
        copy(buf);
        return n_body;
    }
    void writeBinary(FILE* fp) const {
        FileHeader buf;
        buf.copy(*this);
        if ( !fwrite(&buf, sizeof(buf), 1, fp) ) {
            errorMessage("The header has NOT been correctly written.");
            PS::Abort();
        }
    }

    FileHeader(){ n_body = 0; id_next = 0; time = 0.; }
    FileHeader(PS::S32 n_body0,
               PS::S32 id_next0,
               PS::F64 time0,
               Energy e_init0,
               Energy e_now0){
        n_body = n_body0;
        id_next = id_next0;
        time = time0;
        e_init = e_init0;
        e_now = e_now0;
    }
};

#ifdef OUTPUT_DETAIL
template<class Tpsys>
void calcKineticEnergy(const Tpsys & pp,
                       PS::F64 & ekin)
{
    ekin = 0.;   
    PS::F64 ekin_loc = 0.0;    
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    
    for(PS::S32 i = 0; i < n_loc; i++){
        ekin_loc += pp[i].mass * pp[i].vel * pp[i].vel;
    }
    ekin_loc *= 0.5;
    
    ekin = PS::Comm::getSum(ekin_loc);
}
#endif


