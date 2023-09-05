import numpy as np

hoehen = [0,3.,6.,9.]
kraefte = [0,50., 100., 200.]
Ersatzsteifigkeit = 31500000.0


# from suds.client import Client

# client = Client('http://localhost:8081/wsdl')
# new = client.service.get_active_model()+'wsdl'
# model = Client(new)

# print(model.service.get_section(1))

def FEM_ersatzstab(hoehen, kraefte, Ersatzsteifigkeit):
    import os
    import sys
    baseName = os.path.basename(__file__)
    dirName = os.path.dirname(__file__)

    sys.path.append(dirName + r'/../..')
    
    
    from RFEM.enums import NodalSupportType, StaticAnalysisType, NodalLoadDirection
    from RFEM.initModel import Model, Calculate_all
    from RFEM.BasicObjects.material import Material
    from RFEM.BasicObjects.section import Section
    from RFEM.BasicObjects.node import Node
    from RFEM.BasicObjects.member import Member
    from RFEM.TypesForNodes.nodalSupport import NodalSupport
    from RFEM.LoadCasesAndCombinations.staticAnalysisSettings import StaticAnalysisSettings
    from RFEM.LoadCasesAndCombinations.loadCase import LoadCase
    from RFEM.Loads.nodalLoad import NodalLoad
    from RFEM.Results.resultTables import GetMaxValue, GetMinValue, ResultTables
    
    
    # RFEM 6
    Model(True, "Abschaetzung_T1")
    Model.clientModel.service.begin_modification()
    
    
   
    Material(1, 'C30/37')
    
    section_params = {
        'name':'Ersatzquerschnitt',
        'shear_stiffness_deactivated' : True,
        'warping_stiffness_deactivated' : True,
        'thin_walled_model' : True,
        'us_spelling_of_properites' : False,
        'stress_smoothing_to_avoid_singularities' : False,
        'area_axial' : 0.00538,
        'area_shear_y' : 0.002657379815189324,
        'area_shear_z' : 0.0020336126167305626,
        'inclination_principal_axes' : 0.0,
        'rotation_angle' : 0.0,
        'location_of_centroidal_axis_y' : 0.075,
        'location_of_centroidal_axis_z' : 0.15,
        'moment_of_inertia_bending_y' : Ersatzsteifigkeit,
        'moment_of_inertia_bending_z' : 6.04e-06,
        'moment_of_inertia_torsion' : 1.99e-07,
        'depth_temperature_load' : 0.3,
        'width_temperature_load' : 0.15,
        'material' : 1,
        'comment' : None,
        'is_generated' : False,
        'metadata_for_export_import' : None,

        } 
    Section(1,name='Basic', material_no=1, params= section_params)
    
    for hoehe in enumerate(hoehen):
        Node(hoehe[0]+1, hoehe[1], 0.0,0.0)
        
    for i in range(len(hoehen)-1):
        Member(i+1, i+1, i+2, 0.0, 1, 1)
    
    NodalSupport(1, '1', NodalSupportType.FIXED)
    
    StaticAnalysisSettings(
        1, '1. Ordnung', StaticAnalysisType.GEOMETRICALLY_LINEAR)
    
    LoadCase(1, 'Ersatzkraefte',self_weight=[False])
    
    for i in range(len(kraefte)-1):
        NodalLoad(i+1, 1, nodes_no=f'{i+2}', load_direction=NodalLoadDirection.LOAD_DIRECTION_GLOBAL_Z_OR_USER_DEFINED_W, magnitude=kraefte[i+1])


    Calculate_all()
    Deformation_Node_4 = ResultTables.NodesDeformations(loading_no=1, object_no=4)
    
    
    Model.clientModel.service.finish_modification()
    return Deformation_Node_4

FEM_ersatzstab(hoehen=hoehen, kraefte= kraefte, Ersatzsteifigkeit=Ersatzsteifigkeit)

