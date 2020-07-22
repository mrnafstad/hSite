<template>
  <div>
    <img id="profilepic" src="../assets/profilepic.jpg"> <br/>
    <h1>{{ greeting }}!</h1>
    <AboutMe />
  </div>
  
</template>

<script>
import AboutMe from '@/components/AboutMe.vue'
export default {
  components: {
    AboutMe
  },
    data: function() {
      return {
        greeting: "Hello world, it's me"
      };
    }
  };
</script>

<style scoped>
  h1 {
    padding-top: 10px;
    font-size: 2em;
    text-align: center;
    color: red;
  }
  #profilepic {
    border-radius: 50%;
    width: 300px;
    padding: 10px;
}

</style>