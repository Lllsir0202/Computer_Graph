import random
import math

def generate_box(file, id, x, y, z, size):
    """生成一个立方体的6个面"""
    file.write(f"""
Begin Model
Model Box{id}
Translation {x} {y} {z}

Plane BoxPlane1 White
N 0 1 0
P {size} {size} 0
U -{size*2} 0 0
V 0 0 -{size*2}

Plane BoxPlane2 White
N 0 -1 0
P {size} -{size} 0
U -{size*2} 0 0
V 0 0 -{size*2}

Plane BoxPlane3 White
N 1 0 0
P {size} -{size} 0
U 0 {size*2} 0
V 0 0 -{size*2}

Plane BoxPlane4 White
N -1 0 0
P -{size} -{size} 0
U 0 {size*2} 0
V 0 0 -{size*2}

Plane BoxPlane5 White
N 0 0 1
P 0 -{size} {size}
U {size*2} 0 0
V 0 {size*2} 0

Plane BoxPlane6 White
N 0 0 -1
P 0 -{size} -{size}
U {size*2} 0 0
V 0 {size*2} 0

End
""")

def generate_scene(output_path):
    with open(output_path, 'w') as f:
        # 写入材质定义
        f.write("""Begin Material
Material White
Prop diffuseColor RGB 0.725 0.71 0.68
Material Red
Prop diffuseColor RGB 0.63 0.065 0.05
Material Green
Prop diffuseColor RGB 0.14 0.45 0.091
End

Begin Model
Model Wall
Translation 0.0 0.0 2056.0
""")
        
        # 写入康奈尔盒墙面
        walls = """Plane LeftWall Red
N -1.0 0.0 0.0
P 556.0 556.0 556.0
U 0 -1112.0 0
V 0 0 -1112.0

Plane RightWall Green
N 1.0 0.0 0.0
P -556.0 556.0 556.0
U 0 -1112.0 0
V 0 0 -1112.0

Plane TopWall White
N 0.0 -1.0 0.0
P 556.0 556.0 556.0
U -1112.0 0 0
V 0 0 -1112.0

Plane BottomWall White
N 0.0 1.0 0.0
P 556.0 -556.0 556.0
U -1112.0 0 0
V 0 0 -1112.0

Plane BackWall White
N 0.0 0.0 -1.0
P 556.0 556.0 556.0
U -1112.0 0 0
V 0 -1112.0 0

End
"""
        f.write(walls)

        # 生成小立方体阵列
        box_size = 30
        box_id = 0
        for layer in range(4):  # 4层
            base_y = -456 + layer * box_size * 2  # 从底部往上堆积
            for row in range(10):  # 10行
                for col in range(10):  # 10列
                    x = -400 + col * box_size * 2.2  # 间隔排列
                    z = 1400 + row * box_size * 2.2  # 从后向前排列
                    # 添加随机偏移
                    x_offset = random.uniform(-5, 5)
                    z_offset = random.uniform(-5, 5)
                    if layer > 0:  # 上层的方块添加更多随机性
                        x_offset *= 1.5
                        z_offset *= 1.5
                    
                    # 确保方块不会超出康奈尔盒边界
                    x = max(min(x + x_offset, 500), -500)
                    z = max(min(z + z_offset, 1900), 1300)
                    
                    generate_box(f, box_id, x, base_y, z, box_size)
                    box_id += 1

        # 添加光源
        f.write("""
Begin Light
Area TopLight
IRV 95.6768 77.1328 62.1616
P 120 550 2176
U -240 0 0
V 0 0 -240
End
""")

if __name__ == "__main__":
    output_path = "../resource/many_boxes.scn"
    generate_scene(output_path)
    print(f"Generated scene with {4*10*10*6} faces")  # 4层*10行*10列*6面
