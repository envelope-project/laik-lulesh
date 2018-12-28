#include "lulesh.h"

void Domain::re_distribute_data_structures(Laik_Group* new_group, Laik_Partitioning* p_exclusive, Laik_Partitioning* p_halo, Laik_Partitioning* p_overlapping, Laik_Transition *t_to_exclusive, Laik_Transition *t_to_halo, Laik_Transition *t_to_overlapping_init, Laik_Transition *t_to_overlapping_reduce){

#ifdef REPARTITIONING
    m_x.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_y.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_z.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_xd.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_yd.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_zd.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_xdd.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_ydd.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_zdd.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
#endif
    m_fx.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_fy.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_fz.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
    m_nodalMass.migrate(new_group, p_overlapping, p_overlapping, t_to_overlapping_init, t_to_overlapping_reduce);
#ifdef REPARTITIONING
    m_dxx.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_dyy.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_dzz.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
#endif
    m_delv_xi.migrate(new_group, p_exclusive, p_halo, t_to_exclusive, t_to_halo);
    m_delv_eta.migrate(new_group, p_exclusive, p_halo, t_to_exclusive, t_to_halo);
    m_delv_zeta.migrate(new_group, p_exclusive, p_halo, t_to_exclusive, t_to_halo);
#ifdef REPARTITIONING
    m_delx_xi.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_delx_eta.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_delx_zeta.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_e.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_p.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_q.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_ql.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_qq.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_v.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_volo.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    //m_vnew.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_delv.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_vdov.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_arealg.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_ss.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);
    m_elemMass.migrate(new_group, p_exclusive,nullptr, nullptr, nullptr);

#endif
    this -> world = new_group;
}
