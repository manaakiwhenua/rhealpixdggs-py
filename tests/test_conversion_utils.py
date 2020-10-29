"""
This Python 3.3 code tests the ``conversion_utils`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

CHANGELOG:

- David Habgood (DH), 2020-10-27: Initial version.
"""


# Import standard modules.
import unittest
from shapely import wkt

# Import my modules.
from rhealpixdggs.conversion_utils import *


test_geom = """
        MULTIPOLYGON (((148.6675 -35.5725, 148.6675 -35.575, 148.67 -35.575, 148.67 -35.5775, 148.6725 -35.5775,
        148.6725 -35.58, 148.675 -35.58, 148.675 -35.5825, 148.6775 -35.5825, 148.6775 -35.585, 148.685 -35.585,
        148.685 -35.5825, 148.7 -35.5825, 148.7 -35.585, 148.7025 -35.585, 148.7025 -35.5875, 148.705 -35.5875,
        148.705 -35.59, 148.7075 -35.59, 148.7075 -35.595, 148.7075 -35.6, 148.71 -35.6, 148.71 -35.605, 148.7125
        -35.605, 148.7125 -35.6075, 148.715 -35.6075, 148.715 -35.61, 148.7275 -35.61, 148.7375 -35.61, 148.7375
        -35.6025, 148.7425 -35.6025, 148.7425 -35.6, 148.745 -35.6, 148.745 -35.5975, 148.7525 -35.5975, 148.7525
        -35.595, 148.7575 -35.595, 148.7575 -35.5975, 148.7625 -35.5975, 148.7625 -35.6, 148.765 -35.6, 148.765
        -35.6025, 148.7675 -35.6025, 148.7675 -35.605, 148.77 -35.605, 148.77 -35.6075, 148.7725 -35.6075, 148.7725
        -35.61, 148.7775 -35.61, 148.7775 -35.615, 148.78 -35.615, 148.78 -35.6175, 148.78 -35.6225, 148.78 -35.625,
        148.7825 -35.625, 148.7825 -35.63, 148.78 -35.63, 148.78 -35.635, 148.7775 -35.635, 148.7775 -35.6375, 148.775
        -35.6375, 148.775 -35.64, 148.7725 -35.64, 148.7725 -35.6425, 148.77 -35.6425, 148.77 -35.6475, 148.77 -35.65,
        148.7725 -35.65, 148.7725 -35.655, 148.7875 -35.655, 148.7875 -35.6575, 148.79 -35.6575, 148.79 -35.66, 148.7925
        -35.66, 148.7925 -35.6625, 148.7975 -35.6625, 148.7975 -35.665, 148.7975 -35.67, 148.7975 -35.675, 148.795
        -35.675, 148.795 -35.68, 148.7975 -35.68, 148.7975 -35.685, 148.795 -35.685, 148.795 -35.6875, 148.7925
        -35.6875,148.7925 -35.69, 148.7925 -35.695, 148.79 -35.695, 148.79 -35.7025, 148.7925 -35.7025, 148.7925
        -35.7075, 148.795 -35.7075, 148.795 -35.71, 148.8 -35.71, 148.8 -35.7125, 148.8025 -35.7125, 148.8025 -35.715,
        148.8025 -35.7175, 148.8 -35.7175, 148.8 -35.72, 148.7975 -35.72, 148.7975 -35.7225, 148.795 -35.7225, 148.795
        -35.725, 148.79 -35.725, 148.79 -35.7275, 148.7825 -35.7275, 148.78 -35.7275, 148.78 -35.725, 148.7775 -35.725,
        148.7775 -35.7225, 148.775 -35.7225, 148.775 -35.7175, 148.7725 -35.7175, 148.7725 -35.7125, 148.77 -35.7125,
        148.77 -35.71, 148.7675 -35.71, 148.7675 -35.7075, 148.765 -35.7075, 148.765 -35.705, 148.7625 -35.705, 148.7625
        -35.7, 148.76 -35.7, 148.76 -35.6975, 148.7575 -35.6975, 148.7575 -35.695, 148.755 -35.695, 148.755 -35.69,
        148.7525 -35.69, 148.7525 -35.6875, 148.75 -35.6875, 148.75 -35.685, 148.745 -35.685, 148.745 -35.6875, 148.74
        -35.6875, 148.74 -35.685, 148.735 -35.685, 148.735 -35.6825, 148.73 -35.6825, 148.73 -35.68, 148.7275 -35.68,
        148.7275 -35.6725, 148.73 -35.6725, 148.73 -35.6675, 148.7275 -35.6675, 148.7275 -35.665, 148.7175 -35.665,
        148.715 -35.665, 148.715 -35.6675, 148.7125 -35.6675, 148.71 -35.6675, 148.71 -35.665, 148.7075 -35.665,
        148.7075 -35.6625, 148.7025 -35.6625, 148.7025 -35.66, 148.695 -35.66, 148.695 -35.665, 148.6925 -35.665,
        148.6925 -35.6625, 148.69 -35.6625, 148.69 -35.66, 148.6875 -35.66, 148.6875 -35.6575, 148.66 -35.6575,
        148.65 -35.6575, 148.65 -35.655, 148.6475 -35.655, 148.6475 -35.6525, 148.6425 -35.6525, 148.64 -35.6525,
        148.64 -35.65, 148.64 -35.6475, 148.635 -35.6475, 148.635 -35.645, 148.6325 -35.645, 148.6325 -35.6425, 148.63
        -35.6425, 148.63 -35.64, 148.6275 -35.64, 148.625 -35.64, 148.625 -35.6375, 148.6225 -35.6375, 148.6225 -35.635,
        148.62 -35.635, 148.62 -35.6325, 148.62 -35.63, 148.6175 -35.63, 148.6175 -35.625, 148.615 -35.625, 148.615
        -35.6225, 148.6125 -35.6225, 148.6125 -35.62, 148.615 -35.62, 148.615 -35.6175, 148.6175 -35.6175, 148.6175
        -35.6075, 148.6225 -35.6075, 148.6225 -35.61, 148.625 -35.61, 148.625 -35.6125, 148.6275 -35.6125, 148.6275
        -35.61, 148.63 -35.61, 148.63 -35.6075, 148.63 -35.605, 148.6325 -35.605, 148.6325 -35.6, 148.635 -35.6, 148.635
        -35.5975, 148.6375 -35.5975, 148.6375 -35.595, 148.64 -35.595, 148.64 -35.59, 148.64 -35.5875, 148.645 -35.5875,
        148.645 -35.585, 148.6475 -35.585, 148.65 -35.585, 148.65 -35.5825, 148.6525 -35.5825, 148.6525 -35.58, 148.655
        -35.58, 148.655 -35.575, 148.6575 -35.575, 148.6575 -35.5725, 148.665 -35.5725, 148.6675 -35.5725)))"""
test_poly = wkt.loads(test_geom)
# test geom from https://gds.loci.cat/geometry/geofabric2_1_1_contractedcatchment/12104622?_format=text/turtle&_view=geometryview
test_feature = ["Contracted_Catchment_12104622", test_poly]
ground_truth_cells_for_catchment_12104622_at_res_9 = ['R785183285', 'R785183287', 'R785183288',
                                                           'R785183514', 'R785183517', 'R785183518',
                                                           'R785183520', 'R785183521', 'R785183522',
                                                           'R785183523', 'R785183524', 'R785183525',
                                                           'R785183526', 'R785183527', 'R785183528',
                                                           'R785183540', 'R785183541', 'R785183542',
                                                           'R785183543', 'R785183544', 'R785183545',
                                                           'R785183547', 'R785183548', 'R78518355',
                                                           'R785183571', 'R785183572', 'R785183575',
                                                           'R785183580', 'R785183581', 'R785183582',
                                                           'R785183583', 'R785183584', 'R785183585',
                                                      'R785183587', 'R785183588', 'R785184035',
                                                      'R785184037', 'R785184038', 'R785184040',
                                                      'R785184041', 'R785184042', 'R785184043',
                                                      'R785184044', 'R785184045', 'R785184046',
                                                      'R785184047', 'R785184048', 'R785184053',
                                                      'R785184056', 'R785184057', 'R785184058',
                                                      'R78518406', 'R78518407', 'R78518408', 'R785184136',
                                                      'R785184137', 'R785184138', 'R785184146',
                                                      'R78518416', 'R785184170', 'R785184171',
                                                      'R785184173', 'R785184174', 'R785184176',
                                                      'R785184177', 'R785184178', 'R785184277',
                                                      'R785184278', 'R785184283', 'R785184286',
                                                      'R785184287', 'R785184288', 'R7851843', 'R78518440',
                                                      'R785184410', 'R785184411', 'R785184412',
                                                      'R785184413', 'R785184414', 'R785184415',
                                                      'R785184416', 'R785184417', 'R785184418',
                                                      'R785184423', 'R785184424', 'R785184425',
                                                      'R785184426', 'R785184427', 'R785184428',
                                                      'R78518443', 'R78518444', 'R78518445', 'R78518446',
                                                      'R78518447', 'R78518448', 'R785184503', 'R785184504',
                                                      'R785184505', 'R785184506', 'R785184507',
                                                      'R785184508', 'R78518451', 'R78518452', 'R78518453',
                                                      'R78518454', 'R78518455', 'R78518456', 'R78518457',
                                                      'R78518458', 'R785184600', 'R785184601',
                                                      'R785184602', 'R785184604', 'R785184605',
                                                      'R785184610', 'R785184611', 'R785184612',
                                                      'R785184613', 'R785184614', 'R785184615',
                                                      'R785184620', 'R785184621', 'R785184622',
                                                      'R785184623', 'R785184624', 'R785184625',
                                                      'R785184700', 'R785184701', 'R785184702',
                                                      'R785184703', 'R785184704', 'R785184705',
                                                      'R785184707', 'R785184708', 'R785184710',
                                                      'R785184711', 'R785184712', 'R785184713',
                                                      'R785184714', 'R785184715', 'R785184716',
                                                      'R785184717', 'R785184718', 'R78518472',
                                                      'R785184732', 'R785184742', 'R785184750',
                                                      'R785184751', 'R785184752', 'R78518480', 'R78518481',
                                                      'R78518482', 'R785184830', 'R785184831',
                                                      'R785184832', 'R785184834', 'R785184835',
                                                      'R785184837', 'R785184838', 'R78518484', 'R78518485',
                                                      'R785184861', 'R785184862', 'R785184864',
                                                      'R785184865', 'R785184870', 'R785184871',
                                                      'R785184872', 'R785184873', 'R785184874',
                                                      'R785184875', 'R785184876', 'R785184877',
                                                      'R785184880', 'R785184881', 'R785184882',
                                                      'R785184883', 'R785184884', 'R785184885',
                                                      'R785184886', 'R785184887', 'R785184888',
                                                      'R785185300', 'R785185303', 'R785185304',
                                                      'R785185305', 'R785185306', 'R785185307',
                                                      'R785185308', 'R78518533', 'R785185346',
                                                      'R785185360', 'R785185361', 'R785185362',
                                                      'R785185363', 'R785185364', 'R785185366',
                                                      'R785185600', 'R785185603', 'R785185604',
                                                      'R785185606', 'R785185607', 'R785185608',
                                                      'R785185616', 'R785185617', 'R78518563',
                                                      'R785185640', 'R785185641', 'R785185642',
                                                      'R785185643', 'R785185644', 'R785185645',
                                                      'R785185646', 'R785185647', 'R785185648',
                                                      'R785185650', 'R785185653', 'R785185656',
                                                      'R78518566', 'R785185670', 'R785185671',
                                                      'R785185672', 'R785185673', 'R785185674',
                                                      'R785185675', 'R785185676', 'R785185677',
                                                      'R785185678', 'R785185683', 'R785187221',
                                                      'R785187222', 'R785187224', 'R785187225',
                                                      'R785187228', 'R78518800', 'R785188010',
                                                      'R785188011', 'R785188012', 'R785188013',
                                                      'R785188014', 'R785188016', 'R785188017',
                                                      'R785188030', 'R785188031', 'R785188032',
                                                      'R785188034', 'R785188035', 'R785188038',
                                                      'R785188040', 'R785188041', 'R785188042',
                                                      'R785188043', 'R785188044', 'R785188045',
                                                      'R785188046', 'R785188047', 'R785188048',
                                                      'R785188053', 'R785188054', 'R785188056',
                                                      'R785188057', 'R785188062', 'R785188070',
                                                      'R785188071', 'R785188072', 'R785188073',
                                                      'R785188074', 'R785188080']
ground_truth_ordered_compressed_cells_for_catchment_12104622_at_res_9 = ['R7851843', 'R78518352', 'R78518355', 'R78518404',
                                                            'R78518406', 'R78518407', 'R78518408', 'R78518416',
                                                            'R78518440', 'R78518443', 'R78518444', 'R78518445',
                                                            'R78518446', 'R78518447', 'R78518448', 'R78518441',
                                                            'R78518451', 'R78518452', 'R78518453', 'R78518454',
                                                            'R78518455', 'R78518456', 'R78518457', 'R78518458',
                                                            'R78518471', 'R78518472', 'R78518480', 'R78518481',
                                                            'R78518482', 'R78518484', 'R78518485', 'R78518488',
                                                            'R78518533', 'R78518563', 'R78518566', 'R78518564',
                                                            'R78518567', 'R78518800', 'R78518804', 'R785183285',
                                                            'R785183287', 'R785183288', 'R785183514', 'R785183517',
                                                                         'R785183518', 'R785183540', 'R785183541', 'R785183542',
                                                                         'R785183543', 'R785183544', 'R785183545', 'R785183547',
                                                                         'R785183548', 'R785183571', 'R785183572', 'R785183575',
                                                                         'R785183580', 'R785183581', 'R785183582', 'R785183583',
                                                                         'R785183584', 'R785183585', 'R785183587', 'R785183588',
                                                                         'R785184035', 'R785184037', 'R785184038', 'R785184053',
                                                                         'R785184056', 'R785184057', 'R785184058', 'R785184136',
                                                                         'R785184137', 'R785184138', 'R785184146', 'R785184170',
                                                                         'R785184171', 'R785184173', 'R785184174', 'R785184176',
                                                                         'R785184177', 'R785184178', 'R785184277', 'R785184278',
                                                                         'R785184283', 'R785184286', 'R785184287', 'R785184288',
                                                                         'R785184423', 'R785184424', 'R785184425', 'R785184426',
                                                                         'R785184427', 'R785184428', 'R785184503', 'R785184504',
                                                                         'R785184505', 'R785184506', 'R785184507', 'R785184508',
                                                                         'R785184600', 'R785184601', 'R785184602', 'R785184604',
                                                                         'R785184605', 'R785184610', 'R785184611', 'R785184612',
                                                                         'R785184613', 'R785184614', 'R785184615', 'R785184620',
                                                                         'R785184621', 'R785184622', 'R785184623', 'R785184624',
                                                                         'R785184625', 'R785184700', 'R785184701', 'R785184702',
                                                                         'R785184703', 'R785184704', 'R785184705', 'R785184707',
                                                                         'R785184708', 'R785184732', 'R785184742', 'R785184750',
                                                                         'R785184751', 'R785184752', 'R785184830', 'R785184831',
                                                                         'R785184832', 'R785184834', 'R785184835', 'R785184837',
                                                                         'R785184838', 'R785184861', 'R785184862', 'R785184864',
                                                                         'R785184865', 'R785184870', 'R785184871', 'R785184872',
                                                                         'R785184873', 'R785184874', 'R785184875', 'R785184876',
                                                                         'R785184877', 'R785185300', 'R785185303', 'R785185304',
                                                                         'R785185305', 'R785185306', 'R785185307', 'R785185308',
                                                                         'R785185346', 'R785185360', 'R785185361', 'R785185362',
                                                                         'R785185363', 'R785185364', 'R785185366', 'R785185600',
                                                                         'R785185603', 'R785185604', 'R785185606', 'R785185607',
                                                                         'R785185608', 'R785185616', 'R785185617', 'R785185650',
                                                                         'R785185653', 'R785185656', 'R785185683', 'R785187221',
                                                                         'R785187222', 'R785187224', 'R785187225', 'R785187228',
                                                                         'R785188010', 'R785188011', 'R785188012', 'R785188013',
                                                                         'R785188014', 'R785188016', 'R785188017', 'R785188030',
                                                                         'R785188031', 'R785188032', 'R785188034', 'R785188035',
                                                                         'R785188038', 'R785188053', 'R785188054', 'R785188056',
                                                                         'R785188057', 'R785188062', 'R785188070', 'R785188071',
                                                                         'R785188072', 'R785188073', 'R785188074', 'R785188080']


class ConversionUtilsTestCase(unittest.TestCase):
    def test_call_get_finest(self):

        assert(str(call_get_finest(test_poly))) == 'R78518'


    def test_CellZoneFromPoly(self):
        """
        Tests correct cells are output for the Geofabric Contracted Catchment 12104622 at resolution 9, without ordering
        or compression.
        """
        cells_obj_list = CellZoneFromPoly(feature=test_feature, res_limit=9, return_cells=True).cells_list
        self.cell_str_list = [str(cell) for cell in cells_obj_list]
        for cell_str in self.cell_str_list:
            assert cell_str in ground_truth_cells_for_catchment_12104622_at_res_9


    def test_compress_order_cells(self):
        """
        Tests correct cells are output for the Geofabric Contracted Catchment 12104622 at resolution 9, with ordering
        and compression.
        """
        compressed_cells = compress_order_cells(ground_truth_cells_for_catchment_12104622_at_res_9)
        assert compressed_cells == ground_truth_ordered_compressed_cells_for_catchment_12104622_at_res_9

if __name__ == "__main__":
    unittest.main()
