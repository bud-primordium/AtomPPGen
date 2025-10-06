势反演模块 (invert)
====================

.. currentmodule:: atomppgen.invert

概述
----

从 Troullier-Martins 伪轨道反演半局域势 :math:`V_l(r)`。

核心功能
--------

.. autofunction:: invert_semilocal_potential

数据类
------

.. autoclass:: InvertResult
   :members:
   :undoc-members:

内部函数
--------

这些函数通常不需要直接调用，由 ``invert_semilocal_potential`` 内部使用。

.. autofunction:: atomppgen.invert._invert_inner_analytical
.. autofunction:: atomppgen.invert._invert_outer_spline
.. autofunction:: atomppgen.invert._interpolate_nodes
.. autofunction:: atomppgen.invert._count_nodes
