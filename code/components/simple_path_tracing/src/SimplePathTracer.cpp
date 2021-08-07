#include "server/Server.hpp"

#include "SimplePathTracer.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"
#include<stack>

namespace SimplePathTracer
{
    float calDis(Photon a, Vec3 b);
    float vecGetAxis(Vec3 a, int axis);


    RGB SimplePathTracerRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    float random_float() {
        srand(time(NULL));//设置随机数种子，使每次产生的随机序列不同
        return rand() % 100 / (float)100;
    }

    void SimplePathTracerRenderer::generatePhoton(AreaLight a, Vec3& ori, Vec3& dir ) {
        Vec3 normal = glm::cross(a.u, a.v);
        uniform_real_distribution<float> uni;
        ori = Vec3(a.position[0] + random_float() * (a.u[0]), a.position[1], a.position[2] + random_float() * (a.v[2]));
        float mode = 0;
        do {//用循环找一个方向射光出去
            dir = Vec3(2 * uni(eng), 2 * uni(eng), 2 * uni(eng)) - Vec3(1, 1, 1);
            mode = glm::dot(dir, dir);
        } while (mode >= 1.0 || mode == 0.0 || glm::dot(dir, normal) < 0);

    }

    void SimplePathTracerRenderer::generatePhotonMap() {
        Vec3 ori;
        Vec3 dir;
        Vec3 pow(1, 1, 1);
        float powerscale;
        while (phtmap->p_num<phtmap->max) {
            for (auto a : scene.areaLightBuffer) {
                generatePhoton(a, ori, dir);
                //cout << dir<<endl;
                Photon tp(ori, dir, pow);
                tracePhoton(tp, 0);
            }
        }
        //开始balance
        phtmap->balance(phtmap->p_list,0, phtmap->p_num-1,phtmap->kd);
        //nearestNeighbor* np = new nearestNeighbor(Vec3{-35,-20,-64});
        //np->getNearest(phtmap);
        //delete np;
    }

    void SimplePathTracerRenderer::tracePhoton(Photon& p, int depth) {
        if (depth > this->depth)//超过最大搜索深度
            return;
        HitRecord hitObject = closestHitObject(p);
        Ray p_ray(p.origin,p.direction);
        if (hitObject) {
            float cosine = glm::dot(p.direction, hitObject->normal) / sqrt(glm::dot(p.direction, p.direction))/ sqrt(glm::dot(hitObject->normal, hitObject->normal));
            //p.power = p.power * abs(glm::dot(p.direction, hitObject->normal));
            p.power = p.power*abs(cosine);
            Photon tphoton(hitObject->hitPoint, hitObject->normal, p.power);
            //cout << hitObject->hitPoint << endl;
            int rus = russian(0.5,0.8);

            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(p_ray, hitObject->hitPoint, hitObject->normal);

            auto diffuseColor = shaderPrograms[mtlHandle.index()]->mymaterial.getProperty<Property::Wrapper::RGBType>("diffuseColor");
            Vec3 tcolor = (*diffuseColor).value;

            
            switch (rus)
            {
            case diffusion:
                Vec3 temp = tphoton.power;
                tphoton.power = tphoton.power*tcolor;
                tphoton.direction = scattered.ray.direction;
                phtmap->add_p(&tphoton);//碰到漫反射条件，存下来
                if (abs(random_float()) < 0.5) {
                    tphoton.power = temp* scattered.attenuation;
                    tracePhoton(tphoton, depth + 1);
                }
                break;
            case specular:
                tphoton.power = tphoton.power;
                tracePhoton(tphoton,depth+1);
                break;
            default:
                tphoton.power = tphoton.power * tcolor;
                phtmap->add_p(&tphoton);//碰到漫反射条件，存下来
                return;
            }

        }
    }


    void SimplePathTracerRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) { //off 是第几个进程（就是offset），step是总的进程数量（作为步长）
        for(int i=off; i<height; i+=step) {//每个进程负责index是off+i*step的像素
            for (int j=0; j<width; j++) {
                Vec3 color{0, 0, 0};
                Vec3 irra{ 0,0,0 };
                //——————————————————————————————————————————
                 //下面的代码是光子映射算法默认只渲染一次的情况，请根据需求注释或取消注释
                auto r = defaultSamplerInstance<UniformInSquare>().sample2d();// 获得一个2维的随机数
                float rx = r.x; //
                float ry = r.y;//
                float x = (float(j) + rx) / float(width);
                float y = (float(i) + ry) / float(height);
                auto ray = camera.shoot(x, y);
                irra = getIrradiance(ray);//无论采样率多少，只渲染一次
                 //——————————————————————————————————————————


                for (int k=0; k < samples; k++) {
                    auto r = defaultSamplerInstance<UniformInSquare>().sample2d();// 获得一个2维的随机数罢了
                    float rx = r.x; //随机数的第一个值
                    float ry = r.y;// ........第二个值
                    float x = (float(j)+rx)/float(width);
                    float y = (float(i)+ry)/float(height);
                    //相当于上面生成了随机的一个像素点，然后从这个像素点的坐标射出去
                    auto ray = camera.shoot(x, y);

                    //——————————————————————————————————————————
                    // 下面一行代码是simple path tracer，如果需要使用路径追踪，请取消注释
                    //color += trace(ray, 0);//0是递归深度
                    // ——————————————————————————————————————————

                }
                color /= samples;//每次采样结果累加除以采样次数
                color += irra;
                //color = clamp(color);
                color = gamma(color);

                pixels[(height - i - 1) * width + j] = { color, 1 };//这赋值这么奇怪是因为，00在左下角，这样第一个像素赋值就在左上角了
            }
        }
    }

    auto SimplePathTracerRenderer::render() -> RenderResult {
        // shaders
        shaderPrograms.clear();
        ShaderCreator shaderCreator{};
        for (auto& m : scene.materials) {
            shaderPrograms.push_back(shaderCreator.create(m, scene.textures));
        }

        RGBA* pixels = new RGBA[width*height]{};

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);

        this->generatePhotonMap();//构建光子图


        const auto taskNums = 8;/////////////////线程数量
        thread t[taskNums];
        for (int i=0; i < taskNums; i++) {
            t[i] = thread(&SimplePathTracerRenderer::renderTask,
                this, pixels, width, height, i, taskNums);
        }
        for(int i=0; i < taskNums; i++) {
            t[i].join();
        }
        getServer().logger.log("Done...");
        return {pixels, width, height};
    }

    void SimplePathTracerRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord SimplePathTracerRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;
        for (auto& s : scene.sphereBuffer) {
            auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer) {
            auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& p : scene.planeBuffer) {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        return closestHit; 
    }

    HitRecord SimplePathTracerRenderer::closestHitObject(const Photon& p) {
        Ray r;
        r.setDirection(p.direction);
        r.setOrigin(p.origin);
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;
        for (auto& s : scene.sphereBuffer) {
            auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer) {
            auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& p : scene.planeBuffer) {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        return closestHit;
    }


    
    tuple<float, Vec3> SimplePathTracerRenderer::closestHitLight(const Ray& r) {
        Vec3 v = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {});
        for (auto& a : scene.areaLightBuffer) {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t) {
                closest = hitRecord;
                v = a.radiance;
            }
        }
        return { closest->t, v };
    }
    RGB SimplePathTracerRenderer::trace(const Ray& r, int currDepth) {
        if (currDepth == depth) return scene.ambient.constant;//达到递归最大深度，就停止递归
        auto hitObject = closestHitObject(r);
        auto [ t, emitted ] = closestHitLight(r);

        //Vec3 pm = getIrradiance(hitObject, r.direction);
        //Vec3 pm = {0,0,0};

        // hit object
        if (hitObject && hitObject->t < t) { //在打到光源之前打到了物体
                auto mtlHandle = hitObject->material;
                auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
                auto scatteredRay = scattered.ray;
                auto attenuation = scattered.attenuation;
                auto emitted = scattered.emitted;
                auto next = trace(scatteredRay, currDepth + 1);//继续递归
                float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
                float pdf = scattered.pdf;
                /**
                 * emitted      - Le(p, w_0)
                 * next         - Li(p, w_i)
                 * n_dot_in     - cos<n, w_i>
                 * atteunation  - BRDF
                 * pdf          - p(w)
                 **/
                return emitted + attenuation * next * n_dot_in / pdf;//这里是是当前射出的光+上一次衰减射入的光
        }

        // 
        else if (t != FLOAT_INF) {//没打到物体且打到光源
            return emitted;//这样得到的光就是没有衰减的直接从光源得到
        }
        else {
            return Vec3{0};
        }
    }

    RGB SimplePathTracerRenderer::getIrradiance(Ray r) {
        Vec3 ir{ 0,0,0 };
        Vec3 norm = r.direction;
        // hit object
        auto hitObject = closestHitObject(r);

        if (hitObject) { //打到物体了
            nearestNeighbor np(hitObject->hitPoint);


            //------------可选择使用加速或不加速---------------------------
            np.getNearest(phtmap); //使用kdtree
            //--------------------------------------------------------------
            //np.getNeighborwithoutKDtree(phtmap);//不使用kdtree
            //--------------------------------------------------------------


            if (np.found < 2)
                return Vec3{ 0,0,0 };
            
            int p_num = np.found < np.sample ? np.found: np.sample;
            for (int i = 0; i < p_num; i++)
            {
                ir += np.dict[np.diskey[i]].power;
            }
            ir = float(150)*ir / ( 4*PI * np.diskey[np.found - 1]);
            //std::cout << ir<<endl;
            return ir;
        }
        else return Vec3{ 0, 0, 0 };
    }

    void photonMap::add_p(Photon* p) {
        if (p_num < max) {
            p_list[p_num].origin = p->origin;
            p_list[p_num].direction = p->direction;
            p_list[p_num].power = p->power;
            p_num++;
        }
        else cout << "Photon number overflow, not added.\n";
    }// to add photon into my photon list


    void photonMap::balance(Photon* tempPhoton, int start, int end,kdtree* kd) {
        //temphton 是复制下来的用来重组的p_list
        //start是从这次的递归开始的tempphoton的索引
        //end是结束的索引
        //index是这一次要赋值回去给p_list的索引
        if (end==start) {
            kd->val = tempPhoton[start];
            kd->val.axis = -1;
            return;
        }
        //for (int i = start; i <= end; i++)
        //    std::cout << tempPhoton[i].origin  << endl;
        //std::cout << endl;

        Photon* seg = new Photon[end - start + 1];
        for (int i = 0; i < end - start + 1; i++) {
            seg[i] = tempPhoton[i + start];
        }
        int axis =findDim(seg,0,end-start);
        heapSort(seg, 0, end-start, axis);
        //for (int i = 0; i < end - start + 1; i++)
        //    std::cout << seg[i].origin << endl;
        //std::cout << endl;


        for (int i = 0; i < end - start + 1; i++) {
            tempPhoton[i + start] = seg[i];
        }
        //for (int i = 0; i <= p_num; i++)
        //    std::cout<< tempPhoton[i].origin << " axis=" << axis<<endl;
        //std::cout << "__________________________________________________";
        // median
        //std::cout << endl;
        int med = calMedian(start, end); //
        tempPhoton[med].axis = axis; //p_list[med]is the split element and should assign the axis to it at where it splits
        kd->val=tempPhoton[med] ;// assign it to the index it should be;
        //std::cout <<"tree node:"<< kd->val.origin<<endl;
        if (med-1 >= start) {
            //this is left branch 
            kd->left = new kdtree(kd);
            balance(tempPhoton, start, med - 1,kd->left);
        }
        if(med+1<=end) {
            kd->right = new kdtree(kd);
            balance(tempPhoton, med + 1, end,kd->right);
        }
        delete[] seg;
    }

    int photonMap::calMedian(int start, int end)
    {
        int num = end - start + 1;
        return start + num / 2;


        //int as = 1, b = 2;
        //while (as < num)
        //{
        //    as += b;
        //    b *= 2;
        //}
        //if (as == num) {
        //    return start + num / 2;
        //}
        //b /= 2;
        //if (as - b / 2 < num) {
        //    return start + as / 2;
        //}
        //else
        //    return start + as / 2 - (as - b / 2 - num);
    }

    void photonMap::heapSort(Photon* plist, int start, int end,int axis) {
        while (start != end) {
            int checknode = (end - 1) / 2;
            while (true) {
                if (checknode < start)
                    break;
                int leftchild = checknode * 2 + 1;//right =left+1
                float parentval, leftval, rightval;
                switch (axis)
                {
                case 0:
                    parentval = plist[checknode].origin.x;
                    leftval = plist[leftchild].origin.x;
                    rightval = plist[leftchild + 1].origin.x;
                    break;
                case 1:
                    parentval = plist[checknode].origin.y;
                    leftval = plist[leftchild].origin.y;
                    rightval = plist[leftchild + 1].origin.y;
                    break;
                case 2:
                    parentval = plist[checknode].origin.z;
                    leftval = plist[leftchild].origin.z;
                    rightval = plist[leftchild + 1].origin.z;
                    break;

                default:
                    std::cout << "axis is not supposed to be " << axis << endl;
                    break;
                }

                float compare;
                if (leftchild + 1 <= end)
                    compare = std::max(leftval, rightval);
                else compare = leftval;

                if (parentval < compare) { // 如果父亲它的孩子

                    if (leftchild + 1 <= end)//右孩子存在
                    {
                        if (leftval > rightval) { //和左孩子替换
                            //left child is greater
                            swap(plist[leftchild], plist[checknode]);
                        }
                        else {
                            swap(plist[leftchild + 1], plist[checknode]);
                        }
                    }
                    else swap(plist[leftchild], plist[checknode]);//只有左孩子

                    if (leftchild * 2 + 1 <= end)
                        checknode = leftchild;
                }
                else checknode--; //如果大于了，那就看前一个父节点
            }
            std::swap(plist[0],plist[end]);
            //std::cout << plist[end].origin<<endl;
            end--;
        }
    }

    int photonMap::findDim(Photon*plist,int start,int end) { //find the axis that contains in the largest cube
        int axis = 0;
        int dis = 0;
        for (int a=0; a < 3; a++) {
            float maxa = -1000000, mina = 1000000;
            for (int i = start; i <= end; i++) {
                maxa = std::max(plist[i].origin[a],maxa);
                mina = std::min(plist[i].origin[a],mina);
            }
            if (maxa - mina > dis) {
                axis = a;
                dis = maxa - mina;
            }
        }
        return axis;
    }

    void nearestNeighbor::findNeighbor(nearestNeighbor* np, kdtree* kd) {
        if (kd->left)
        {
            float dis = vecGetAxis(np->pos, kd->val.axis) - vecGetAxis(kd->val.origin, kd->val.axis);
            if (dis < 0) {
                findNeighbor(np, kd->left);
                if (dis*dis < np->l2dis && kd->right!=nullptr)
                    findNeighbor(np, kd->right);
            }
            else if(kd->right) {
                findNeighbor(np, kd->right);
                if (dis*dis < np->l2dis && kd->left!=nullptr)
                    findNeighbor(np, kd->left);
            }
        }
        float dis2 = glm::dot(np->pos - kd->val.origin, np->pos - kd->val.origin);

        //std::cout << dis2<<endl;
        if (dis2 > np->l2dis)
            return;
        else if(found<max_photons) {
            dict[dis2] = kd->val;
            diskey[found++] = dis2;
        }
        //else {
        //    if (np->isfull) {
        //        for (int i = np->found >> 1; i >= 1; i--) {
        //            int par = 1;
        //            Photon* tmp = np->photons[i];
        //            float* tmp_dist2 = np->l2dis[i];
        //        }
        //    }
        //}

    }


    float calDis(Photon a, Vec3 b) {
        Vec3 t = a.origin - b;
        return glm::dot(t, t);
    }
    float vecGetAxis(Vec3 a, int axis)
    {
        switch (axis) {
        case 0:
            return a.x;
        case 1:
            return a.y; 
        case 2:
            return a.y;
        default:
            return 0; 
        }

    }
    void nearestNeighbor::getNearest(photonMap *phtmap) {
        this->findNeighbor(this, phtmap->kd);
        sort(this->diskey, this->diskey + this->found);
        //for (int i = 0; i < this->sample; i++) {
        //    this->photons[i] = this->dict[this->diskey[i]];
        //    //cout << this->photons[i].origin<<endl;
        //}
    }
    

    int SimplePathTracerRenderer::russian(float diff, float spec) {
        float t;
        srand(time(NULL));//设置随机数种子，使每次产生的随机序列不同
        t = rand() % 100 / (float)100;
        if (t <= diff) // 0<t<dif
        {
            return diffusion;
        }
        else if (t <= spec) {//dif<t<spec
            return specular;
        }
        else {//spec<i<1
            return absorb;
        }

    }


    void nearestNeighbor::getNeighborwithoutKDtree(photonMap* phtmap) {
        phtmap->p_list;
        for (int i = 0; i < phtmap->max; i++)
        {
            if (this->found >= this->max_photons)
                break;
            Vec3 ppos = phtmap->p_list[i].origin-this->pos;
            float dis = glm::dot(ppos, ppos);
            if (dis < this->l2dis)
            {
                this->dict[dis] = phtmap->p_list[i];
                this->diskey[found] = dis;
                this->found++;
            }
        }
        sort(this->diskey, this->diskey + this->found);
    }


    //float SimplePathTracerRenderer::getIrradiance(nearestNeighbor* n)
    //{
    //    if (n->found < n->sample)
    //        return;
    //    n->getNearest(phtmap);
    //    float r = n->diskey[n->sample];//sample 是最远的点的index；
    //    Vec3 brdf{ 1,1,1 };
    //    float flux = 0;
    //   
    //    for (int i = 0; i < n->sample; i++)
    //        flux += glm::dot(brdf, n->dict[n->diskey[i]].origin);
    //    float lr = flux / (2 * PI * r * r);
    //    return lr;
    //}
}






















    //std::stack<kdtree*> path;
    //kdtree* troot = this;

    //// 先找到target所属的子空间，得到path+----------------------------------------------------
    //while (this->left != nullptr && this->right != nullptr) {
    //    path.push(troot);
    //    int axis = troot->val.axis;
    //    float split;
    //    float tempnum;// target 对应维度上的值
    //    switch (axis) {
    //    case 0:
    //        split = troot->val.origin.x;
    //        tempnum = target.x;
    //        break;
    //    case 1:
    //        split = troot->val.origin.y;
    //        tempnum = target.y;
    //        break;
    //    case 2:
    //        split = troot->val.origin.z;
    //        break;
    //    default:
    //        split = troot->val.origin.z;
    //        tempnum = target.z;
    //        break;
    //    }
    //    troot = tempnum <= split ? troot->left : troot->right;// 大于则往左子树去，小于则往左子树去
    //}
    //path.push(troot);
    ////------------------------------------------------------------------------------------------------
    //kdtree* tparent;
    //// 开始回溯并且找n个点
    //for (kdtree* k = path.top(); !path.empty(); k = path.top()) {
    //    tparent = k->parent;
    //    float dis = calDis(k->val, target);
    //    float split = tparent->getAxis(tparent->val.axis);
    //    if (abs(vecGetAxis(target, k->val.axis) - split) < dis)



    //}

    //}
