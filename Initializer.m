classdef Initializer < handle
    properties
        dvec DesignVector
        wingDesign WingDesign
        mda MDA
        optimizer Optimizer
        CD_wing_initial
        CL_initial
        q_design_initial
        drag_fus_initial

        V_MO_initial 
        
        wing_tank_volume_initial
        internal_tank_volume
        S_initial
        c_root_initial
        c_kink_initial
        c_tip_initial

    end
    methods
        function obj = Initializer(dvec)
            arguments
                dvec DesignVector
            end
            obj.dvec = dvec;
            obj.wingDesign = WingDesign(obj.dvec);
            obj.mda = MDA(obj.wingDesign,Const.W_TO_max_initial,Const.W_ZF_initial);
            obj.optimizer = Optimizer(obj.dvec);
            
            [obj.CL_initial, obj.CD_wing_initial] = obj.optimizer.calcCL_CD(Const.W_TO_max_initial,obj.optimizer.wingDesign.W_fuel);
            CDaminusw = obj.CL_initial/Const.LD_initial - obj.CD_wing_initial;
            q_design_initial = obj.optimizer.calculateDesignDynamicPressure();
            obj.drag_fus_initial = CDaminusw * q_design_initial * obj.wingDesign.S*2;
            obj.V_MO_initial = 0.86 * obj.wingDesign.a;
            obj.wing_tank_volume_initial = obj.wingDesign.calculateWingTankVolume();
            obj.internal_tank_volume = max(0,Const.fuel_weight_max_ref/(0.81715e3)-obj.wing_tank_volume_initial);
            obj.S_initial = obj.wingDesign.S;
            obj.c_root_initial = obj.wingDesign.c_root;
            obj.c_kink_initial = obj.wingDesign.c_kink;
            obj.c_tip_initial = obj.wingDesign.c_tip;
            
        end
    end
end