from carousel_conf import conf
import rawe

if __name__=='__main__':
    autogenDir = 'autogen'
    topname = 'kite'
    dae = rawe.models.crosswind_drag_mode(conf)
    rawe.utils.mkprotos.writeAll(dae, topname, autogenDir,haskellDirs=['plot-ho-matic'])

